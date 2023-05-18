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

    pos_junc_f = '../src/2_GET_REF_JUNCS_REFSEQ/BAM_REF_Intersection/'+str(SEQ_LEN)+"bp/100_juncs/d_a_paralog_removed.bed"
    pos_MANE_f = "../src/5_GET_REF_JUNCS_MANE/BAM_REF_Intersection/"+str(SEQ_LEN)+"bp/100_juncs/d_a_paralog_removed.bed"
    pos_ALTS_f = "../src/6_GET_REF_JUNCS_REFSEQ_ALTS/BAM_REF_Intersection/"+str(SEQ_LEN)+"bp/100_juncs/d_a_paralog_removed.bed"

    neg_1_junc_f = '../src/4_GET_BAM_NEG_OPP_STRAND_JUNCS/BAM_junctions/'+str(SEQ_LEN)+"bp/1_juncs/d_a.bed"
    neg_random_junc_f = '../src/3_GET_BAM_NEG_OPP_STRAND_JUNCS_RANDOM/NEG_rev_junctions/'+str(SEQ_LEN)+"bp/d_a/d_a.bed"

    junc_fs = [pos_junc_f, pos_MANE_f, pos_ALTS_f, neg_1_junc_f, neg_random_junc_f]

    output_dir = "./dataset/"
    output_files = [output_dir+"pos/", output_dir+"pos_MANE/", output_dir+"pos_ALTS/", output_dir+"neg_1/", output_dir+"neg_random/"]

    for output_file in output_files:
        os.makedirs(output_file, exist_ok=True)

    nums = [10000, 10000, 10000, 10000, 10000]

    COUNTER = 0
    global_df = pd.DataFrame(columns = [0, 1, 2, 3, 4, 5])

    for junc_fidx in range(0, len(junc_fs), 1):
        print("junc_fidx: ", junc_fidx)
        junc_f = junc_fs[junc_fidx]
        print("junc_f          : ", junc_f)

        os.makedirs(output_files[junc_fidx]+"splam/", exist_ok=True)
        os.makedirs(output_files[junc_fidx]+"spliceai/", exist_ok=True)

        junc_df = pd.read_csv(junc_f, delimiter="\t", header=None)
        # Selecting junctions only on chr1 and chr9 (testing dataset).
        junc_df = junc_df.loc[((junc_df[0] == "chr1") | (junc_df[0] == "chr9"))]
        junc_df = junc_df.loc[junc_df[1] > 0]

        print("nums[junc_fidx]: ", nums[junc_fidx])
        print("len(junc_df)   : ", len(junc_df))
        if nums[junc_fidx] <= len(junc_df):
            junc_df = junc_df.sample(n=nums[junc_fidx], random_state=1).reset_index(drop=True)
        else:
            junc_df = junc_df.sample(n=len(junc_df), random_state=1).reset_index(drop=True)
            # junc_df = junc_df.sample(n=nums[junc_fidx], random_state=1, replace=True).reset_index(drop=True)


        if junc_fidx == 0:
            junc_df[6] = 1
        else:
            junc_df[6] = 0
        print("SPLAM junc_df   : ", junc_df)
        global_df = junc_df
        global_df[2] -= 1

        ################################
        # SPLAM test data curatiolsn
        ################################
        global_df.to_csv(output_files[junc_fidx]+"splam/splam.juncs.bed", sep="\t", header=None, index=0)

        ################################
        # SpliceAI test data curatiolsn
        ################################
        global_df_spliceai = global_df.copy()
        global_df_spliceai[1] -= 5200
        global_df_spliceai[2] += 5200

        print("SpliceAI junc_df: ", global_df_spliceai)
        global_df_spliceai.to_csv(output_files[junc_fidx]+"spliceai/spliceai.juncs.bed", sep="\t", header=None, index=0)

if __name__ == "__main__":
    main()