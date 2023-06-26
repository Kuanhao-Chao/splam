import pandas as pd
import os

def main():
    SEQ_LEN = "800"
    bam_juncs = pd.read_csv('../ALL_JUNCS/junctions.bed', sep="\t", header=None)

    print("Printing all juncs from the BAM file")
    print(bam_juncs)

    neg_1_juncs = pd.read_csv('./NEG_rev_junctions/'+SEQ_LEN+'bp/d_a/d_a.bed', sep="\t", header=None)
    print("Printing all juncs from the neg_1 samples")
    print(neg_1_juncs)

    intersect_df = pd.merge(neg_1_juncs, bam_juncs, how ='inner', on =[0, 1, 2, 5])
    print("intersect_df: ", intersect_df)

if __name__ == "__main__":
    print("Hello world")
    main()