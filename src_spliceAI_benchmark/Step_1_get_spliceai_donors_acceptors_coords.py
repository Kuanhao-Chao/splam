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

# SAMPLE_NUM = 1261186
# SEQ_LENGTH="800"
# QUATER_SEQ_LEN = int(SEQ_LENGTH) // 4

# hg38_ref = "../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa"
# output_bed = "../NEG_junctions/"+SEQ_LENGTH+"bp/neg_junctions.bed"
# output_file = "../INPUTS/"+SEQ_LENGTH+"bp/input_neg.fa"

def main():
    SEQ_LEN = 800
    THRESHOLD = 100

    pos_junc_f = '../src/BAM_REF_Intersection/'+str(SEQ_LEN)+"bp/100_juncs/d_a.bed"
    neg_can_junc_f = '../src/NEG_junctions/'+str(SEQ_LEN)+"bp/d_a.bed"
    neg_noncan_junc_f = '../src/NEG_noncan_junctions/'+str(SEQ_LEN)+"bp/d_a.bed"
    neg_1_junc_f = '../src/BAM_junctions/'+str(SEQ_LEN)+"bp/1_juncs/d_a.bed"
    junc_fs = [pos_junc_f, neg_can_junc_f, neg_noncan_junc_f, neg_1_junc_f]

    output_file = "./OUTPUT/"

    output_files = ["./OUTPUT/pos/", "./OUTPUT/neg_can/", "./OUTPUT/neg_noncan/", "./OUTPUT/neg_1/"]
    for output_file in output_files:
        os.makedirs(output_file, exist_ok=True)

    nums = [3000, 1000, 1000, 1000]

    COUNTER = 0
    global_df = pd.DataFrame(columns = [0, 1, 2, 3, 4, 5])

    for junc_fidx in range(0, 4, 1):
        junc_f = junc_fs[junc_fidx]
        # print(junc_f)
        # print("junc_f: ", junc_f)
        print("output_file : ", output_file)

        junc_df = pd.read_csv(junc_f, delimiter="\t", header=None)

        # junc_df = junc_df.drop(columns=[6, 7])

        junc_df = junc_df.loc[((junc_df[0] == "chr1") | (junc_df[0] == "chr9"))]

        junc_df = junc_df.loc[junc_df[1] > 0]

        junc_df = junc_df.sample(n=nums[junc_fidx]).reset_index(drop=True)

        if junc_fidx == 0:
            junc_df[6] = 1
        else:
            junc_df[6] = 0
        print("junc_df: ", junc_df)

        global_df = junc_df
        # pd.concat([global_df, junc_df], axis=0)
        global_df[2] -= 1
        global_df.to_csv(output_files[junc_fidx]+"splam.juncs.bed", sep="\t", header=None, index=0)

        ################################
        # SpliceAI test data curation
        ################################
        global_df_spliceai = global_df.copy()

        global_df_spliceai[1] -= 5200
        global_df_spliceai[2] += 5200

        print(global_df_spliceai)
        global_df_spliceai.to_csv(output_files[junc_fidx]+"spliceai.juncs.bed", sep="\t", header=None, index=0)

if __name__ == "__main__":
    main()