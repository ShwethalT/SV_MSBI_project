import pandas as pd
import numpy as np
import re

def chrom_downsample(df: pd.DataFrame,chrom_col="CHROM",label_col="label",
major_label=1, minor_label=0, ratio_maj_min=1.0, random_state=42):
    df = df.copy()

    # Split classes
    maj = df[df[label_col] == major_label]
    mino = df[df[label_col] == minor_label]

    # Current counts
    n_min = len(mino)

    # If already small enough, just return
    if len(maj) <= n_min:
        return pd.concat([maj, mino], ignore_index=True).sample(frac=1, random_state=random_state)

    # Compute a global keep-fraction for majority
    keep_frac = round(n_min / len(maj),2)

    # Apply the SAME keep-fraction within each chromosome → preserves CHROM distribution
    pieces = []
    for chrom, g in maj.groupby(chrom_col):
        k = int(np.floor(keep_frac * len(g)))
        if k > 0:
            pieces.append(g.sample(n=k, random_state=random_state))
    maj_down = pd.concat(pieces, ignore_index=True) if pieces else maj.iloc[0:0]

    out = pd.concat([maj_down, mino], ignore_index=True).sample(frac=1, random_state=random_state)
    return out

def ref_with_alt(seq, alt, num=50):
    seq = seq.upper()
    return seq[:num] + alt.upper() + seq[num+1:]

  
def count_homo(seq):
    seq = seq.upper()          
    max_num = 1                
    start = 1                    
    for i in range(1, len(seq)):         
        if seq[i] == seq[i-1]:           
            start += 1
            max_num = max(max_num, start)
        else:                             
            start = 1
    return max_num 