import pandas as pd
import os

def main():
    SEQ_LEN = "800"

    # bam_juncs = pd.read_csv('../1_GET_BAM_JUNCS/BAM_junctions/junctions_100.bed', 
    bam_juncs = pd.read_csv('../ALL_JUNCS/junctions.bed', sep="\t", header=None)

    print("Printing all juncs from the BAM file")
    print(bam_juncs)

    neg_1_juncs = pd.read_csv('./NEG_rev_junctions/'+SEQ_LEN+'bp/d_a/d_a.bed', sep="\t", header=None)
    print("Printing all juncs from the neg_1 samples")
    print(neg_1_juncs)

    intersect_df = pd.merge(neg_1_juncs, bam_juncs, how ='inner', on =[0, 1, 2, 5])
    print("intersect_df: ", intersect_df)
    # out_df = intersect_df[[0, 1, 2, "3_x", "4_x", 5]]
    # ref_juncs = ref_juncs.rename(columns={0:"chr",1:"start", 2:"end", 3:"junc", 4:"score", 5:"strand"})
    # print("out_df: ", out_df)

    
    # # Calling merge() function
    # os.makedirs('./NEG_rev_junctions/'+SEQ_LEN+'bp/', exist_ok=True)
    # d_a_out = './NEG_rev_junctions/'+SEQ_LEN+'bp/d_a.bed'
    # d_out = './NEG_rev_junctions/'+SEQ_LEN+'bp/donor.bed'
    # a_out = './NEG_rev_junctions/'+SEQ_LEN+'bp/acceptor.bed'

    # # intersect_df = pd.merge(ref_juncs, bam_juncs, how ='inner', on =[0, 1, 2, 5])
    # # print("intersect_df: ", intersect_df)
    # # out_df = intersect_df[[0, 1, 2, "3_x", "4_x", 5]]

if __name__ == "__main__":
    print("Hello world")
    main()