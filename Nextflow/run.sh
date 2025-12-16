#!/bin/bash
#SBATCH --job-name=nf_pipeline
#SBATCH --output=nextflow.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=8G

#export PATH="$HOME/miniconda3/bin:$PATH"
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate nf

cd /scratch/home/strikannad/pipeline


nextflow run main.nf -resume --datatype fastq --tumor "$2" --normal "$1"


conda deactivate
