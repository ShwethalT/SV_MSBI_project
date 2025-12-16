nextflow.enable.dsl=2

// process fastq files , output files under directory nanoplot_samplename
process nanoplot {
    publishDir "$projectDir/output/qc", mode: 'copy'
    conda "$projectDir/envs/nanopack.yaml"

    input:
    tuple val(role), val(sample), path(fq)

    output:
    tuple val(role), val(sample), path("nanoplot_${sample}", type: 'dir')

    cpus params.threads

    """
    NanoPlot --fastq ${fq} -o nanoplot_${sample} --threads ${task.cpus} 
    """
 }

//check for adapters and remove them, outputs files as samplename_porechop.fastq.gz

process porechop {
    label 'big_mem'
    publishDir "$projectDir/output/cho", mode: 'copy'
    conda "$projectDir/envs/nanopack.yaml"
   
    input:
    tuple val(role), val(sample) , path(fq) 

    output:
    tuple val(role), val(sample), path("${sample}_porechop.fastq.gz")

    cpus params.threads
	
    """ 
    porechop -i ${fq} -o ${sample}_porechop.fastq.gz --threads ${task.cpus}
    """
}

//trims input fastq and removes poor quality and reads below a threhsold length,outputs file as samplename_chopper.fastq.gz	
process chopper {
    label 'big_mem'
    publishDir "$projectDir/output/cho", mode: 'copy'
    conda "$projectDir/envs/nanopack.yaml"

    input:
    tuple val(role), val(sample) , path(fq)

    output:
    tuple val(role), val(sample), path("${sample}_chopper.fastq.gz")

    cpus params.threads

    """
    chopper -q 7 -l 500 --threads ${task.cpus} -i ${fq} > chopper_out.fastq

    mv chopper_out.fastq "${sample}_chopper.fastq" 

    gzip -f ${sample}_chopper.fastq
    """
}

//qc of filtered fastq, outputs under directory nanoplot_samplename_filt
process filt_nanoplot {
    publishDir "$projectDir/output/qc", mode: 'copy'
    conda "$projectDir/envs/nanopack.yaml"
      
    input:
    tuple val(role), val(sample), path(fq)

    output:
    tuple val(role), val(sample), path("nanoplot_${sample}_filt" , type: 'dir')
    
    cpus params.threads

    """
    NanoPlot --fastq ${fq} -o "nanoplot_${sample}_filt" --threads ${task.cpus}
    """
}

//mapping of fastq with reference file + coord sorting and indexing bams, outputs file as samplename.bam, samplename.bam.bai
process mapping {
    label 'big_mem'
    publishDir "$projectDir/output/mapping", mode: 'copy'
    cache 'lenient'
   
    conda "$projectDir/envs/minimap2_sam.yaml"

    input:
    tuple val(role), val(sample) , path(fq)

    output:
    tuple val(role), val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
       
    cpus params.threads

    """
    minimap2 -ax map-ont -t ${task.cpus} ${params.reference} ${fq} | samtools sort -@ ${task.cpus} -O bam -o ${sample}.bam
    samtools index -@ ${task.cpus} ${sample}.bam
    """
} 

//process checks mapping stats outputs a text file samplename_flagstat.txt
process flagstat {	
    publishDir "$projectDir/output/flagstat", mode: 'copy'
    conda "$projectDir/envs/samtools.yaml"

    input:
    tuple val(role), val(sample) , path(bam) , path(bai)

    output:
    tuple val(role), val(sample), path("${bam.baseName}_flagstat.txt")

    cpus params.threads

    """
    samtools flagstat -@ ${task.cpus} ${bam} > flagstat.txt
    mv flagstat.txt ${bam.baseName}_flagstat.txt
    """
}

//somatic variant calling of tumor + normal bams, outputs samplename.vcf.gz
process somatic {
    label 'big_mem'
    publishDir "$projectDir/output/somatic", mode: 'copy'

    input:
    tuple val(role_t), val(sample_t), path(t_bam), path(t_bai) 
    tuple val(role_n), val(sample_n), path(n_bam), path(n_bai) 

    output:
    path "${sample_t}_${sample_n}/output.vcf.gz"
    path "${sample_t}_${sample_n}/output.vcf.gz.tbi"

    cpus params.threads

    """
    conda run -n clairs /scratch/home/strikannad/pipeline/ClairS/run_clairs -T ${t_bam} -N ${n_bam} -R ${params.reference} -o ${sample_t}_${sample_n} -t ${task.cpus} -p ont_r10_dorado_sup_5khz_ssrs --snv_min_af 0.01 --snv_min_qual 15 --min_coverage 1 --bed_fn ${params.bed}
    """	
}

workflow {
    Channel.fromPath(params.tumor)
	.map{ f -> tuple('T', f.baseName.replaceFirst(/\.f(ast)?q(\.gz)?$/, ''), f) }
	.set {t_fq}
    Channel.fromPath(params.normal)
	.map { f -> tuple('N', f.baseName.replaceFirst(/\.f(ast)?q(\.gz)?$/, ''), f) }
	.set {n_fq}
    first_all_fq = t_fq.mix(n_fq)

    nanoplot(first_all_fq)

    pore = porechop(n_fq)
    
    all_fq = t_fq.mix(pore)

    chop = chopper(all_fq)
    
    filt_nanoplot(chop)
    
    map = mapping(chop)
    
    flagst = flagstat(map)

    t_bam = map.filter { role, sample, bam, bai -> role == 'T' }
    n_bam = map.filter { role, sample, bam, bai -> role == 'N' }

    somatic(t_bam, n_bam)
}
