from Bio import SeqIO
import random
import os
import sys
from pathlib import Path
import pandas as pd

def get_hg38_chrom_size():
    f_chrs = open("../src/hg38_chrom_size.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = int(eles[1])
    return chrs

chrs = get_hg38_chrom_size()

def main():
    SEQ_LEN = 800
    THRESHOLD = 100


    pos_junc_f = '../src/2_GET_REF_JUNCS/BAM_REF_Intersection/'+str(SEQ_LEN)+"bp/100_juncs/d_a.bed"
    pos_junc_f = '../src/2_GET_REF_JUNCS/BAM_REF_Intersection/'+str(SEQ_LEN)+"bp/100_juncs/d_a.bed"
    
    # neg_can_junc_f = '../src/3_GET_CANONICAL_NEG_JUNCS/NEG_can_junctions/'+str(SEQ_LEN)+"bp/d_a.bed"
    # neg_noncan_junc_f = '../src/4_GET_NONCANONICAL_NEG_JUNCS/NEG_noncan_junctions/'+str(SEQ_LEN)+"bp/d_a.bed"
    neg_1_junc_f = '../src/6_GET_BAM_NEG_OPP_STRAND_JUNCS/BAM_junctions/'+str(SEQ_LEN)+"bp/1_juncs/d_a.bed"
    neg_1_junc_random_f = "../src/5_GET_BAM_NEG_OPP_STRAND_JUNCS_RANDOM/NEG_rev_junctions/"+str(SEQ_LEN)+"bp/d_a/d_a.bed"

    # neg_1_junc_f = '../src/6_GET_BAM_NEG_OPP_STRAND_JUNCS/BAM_junctions/'+str(SEQ_LEN)+"bp/1_juncs/d_a.bed"

    # junc_f = '/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/TEST/junction.bed'

    # junc_fs = [pos_junc_f, neg_can_junc_f, neg_noncan_junc_f, neg_1_junc_f]
    # junc_fs = [pos_junc_f, neg_can_junc_f, neg_noncan_junc_f, neg_1_junc_f]

    junc_fs = [pos_junc_f, neg_1_junc_f, neg_1_junc_random_f]

    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"neg_1/", output_dir+"neg_1_random/"]


    for output_file in output_files:
        os.makedirs(output_file, exist_ok=True)

    # Testing [5000, 2500, 2750] samples
    # nums = [1000, 500, 550, 1000]
    # nums = [10000, 1000, 1000, 10000]
    nums = [10000, 10000, 10000]

    COUNTER = 0
    global_df = pd.DataFrame(columns = [0, 1, 2, 3, 4, 5])

    for junc_fidx in range(0, len(junc_fs), 1):
        print("junc_fidx: ", junc_fidx)
        junc_f = junc_fs[junc_fidx]
        print("output_file     : ", output_file)
        print("junc_f          : ", junc_f)
        print("nums[junc_fidx] : ", nums[junc_fidx])

        os.makedirs(output_files[junc_fidx]+"splam/", exist_ok=True)
        os.makedirs(output_files[junc_fidx]+"spliceai/", exist_ok=True)

        junc_df = pd.read_csv(junc_f, delimiter="\t", header=None)
        # Selecting junctions only on chr1 and chr9 (testing dataset).
        junc_df = junc_df.loc[((junc_df[0] == "chr1") | (junc_df[0] == "chr9"))]
        junc_df = junc_df.loc[junc_df[1] > 0]
        junc_df = junc_df.sample(n=nums[junc_fidx]).reset_index(drop=True)

        if junc_fidx == 0:
            junc_df[6] = 1
        else:
            junc_df[6] = 0
        print("junc_df: ", junc_df)

        global_df = junc_df
        global_df[2] -= 1
        global_df.to_csv(output_files[junc_fidx]+"splam/splam.juncs.bed", sep="\t", header=None, index=0)

        ################################
        # SpliceAI test data curation
        ################################
        global_df_spliceai = global_df.copy()

        global_df_spliceai[1] -= 5200
        global_df_spliceai[2] += 5200

        print(global_df_spliceai)
        global_df_spliceai.to_csv(output_files[junc_fidx]+"spliceai/spliceai.juncs.bed", sep="\t", header=None, index=0)

if __name__ == "__main__":
    main()