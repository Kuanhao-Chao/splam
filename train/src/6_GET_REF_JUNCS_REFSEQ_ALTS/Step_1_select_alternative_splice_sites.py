import os 
import re
import pandas as pd

def chr_name_convert():
    f_chrs = open("../../Dataset/Refseq_2_UCSU_chromosome_names.tsv", "r")
    lines = f_chrs.read().splitlines()
    chrs = {}
    for line in lines:
        eles = line.split("\t")
        chrs[eles[0]] = eles[1]
    return chrs

def main():
    JUNC_COUNTER = 0
    THRESHOLD = "100"
    SEQ_LEN="800"
    QUATER_SEQ_LEN = int(SEQ_LEN) // 4

    d_a_out = "./REF_junctions/ref_d_a.bed"
    os.makedirs("./REF_junctions/", exist_ok=True)
    # fw = open("./REF_junctions/ref_d_a.bed", 'w')
    # chrs = chr_name_convert()

    MANE_juncs = pd.read_csv('../5_GET_REF_JUNCS_MANE/BAM_REF_Intersection/800bp/100_juncs/d_a.bed', sep="\t", header=None)

    refseq_juncs = pd.read_csv('../2_GET_REF_JUNCS_REFSEQ/BAM_REF_Intersection/800bp/100_juncs/d_a.bed', sep="\t", header=None)

    print(MANE_juncs)
    print(refseq_juncs)

    intersect_df = pd.merge(refseq_juncs, MANE_juncs, how ='left', on =[0, 1, 2, 5])

    intersect_df = intersect_df[intersect_df.isnull().any(axis=1)]

    print(intersect_df)

    # intersect_df=intersect_df.dropna()
    intersect_df = intersect_df.drop(['3_y', '4_y'], axis=1)
    out_df = intersect_df.rename(columns={0:"chr",1:"start", 2:"end", "3_x":"junc", "4_x":"score", 5:"strand"})
    out_df.to_csv(d_a_out, sep="\t", index=False, header=None)


if __name__ == "__main__":
    main()