import pandas as pd

def main():
    bam_juncs = pd.read_csv('./BAM_junctions/d_a_all_5.bed', sep="\t", header=None)
    neg_juncs = pd.read_csv('./NEG_junctions/pre_filter/d_a.bed', sep="\t", header=None)
    print(bam_juncs)
    print(neg_juncs)
    # Calling merge() function

    d_a_out = './NEG_junctions/d_a.bed'
    d_out = './NEG_junctions/donor.bed'
    a_out = './NEG_junctions/acceptor.bed'

    int_df = pd.merge(neg_juncs, bam_juncs, how ='left', on =[0, 1, 2, 5])
    out_df = int_df[~pd.notna(int_df["3_y"])][[0, 1, 2, "3_x", "4_x", 5]]
    out_df = out_df.rename(columns={0:"chr",1:"start", 2:"end", "3_x":"junc", "4_x":"score", 5:"strand"})
    print(out_df)
    donor_df = out_df.copy()
    acceptor_df = out_df.copy()

    print(donor_df)
    donor_df["start"] = donor_df[["start"]] - 200
    donor_df["end"] = donor_df[["start"]] + 400
    print(donor_df)

    print(acceptor_df)
    acceptor_df["start"] = acceptor_df[["end"]] - 200 -1
    acceptor_df["end"] = acceptor_df[["start"]] + 400
    acceptor_df["junc"] = "JUNC_acceptor"
    print(acceptor_df)
    out_df.to_csv(d_a_out, sep="\t", index=False, header=None)
    donor_df.to_csv(d_out, sep="\t", index=False, header=None)
    acceptor_df.to_csv(a_out, sep="\t", index=False, header=None)

    # print(int_df[[0, 1, 2, "3_x", "4_x", 5]])

    # neg_juncs = neg_juncs.isin(int_df)
    # print("neg_juncs: ", neg_juncs)
    # print(neg_juncs[neg_juncs.isin(int_df)])
    # # print("neg_juncs[1,2]: ", neg_juncs[[1,2]])
    # cond = neg_juncs[[0,1,2]].isin(int_df[[0,1,2]])
    # print(cond)
    # neg_juncs.drop(neg_juncs[cond].index, inplace = True)
    # print(neg_juncs)

if __name__ == "__main__":
    main()