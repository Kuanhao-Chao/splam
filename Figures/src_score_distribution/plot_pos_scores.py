import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import os

output_file = "neg_1"


TOOL = "SPLICEAI"
# TOOL = "SPLAM"

if TOOL == "SPLAM":
    TYPE = "noshuffle"
    d_score_tsv_f = "../../src_tools_evaluation/splam_result/"+output_file+"/splam_all_seq.score.d."+TYPE+"."+output_file+".tsv"
    a_score_tsv_f = "../../src_tools_evaluation/splam_result/"+output_file+"/splam_all_seq.score.a."+TYPE+"."+output_file+".tsv"
    n_score_tsv_f = "../../src_tools_evaluation/splam_result/"+output_file+"/splam_all_seq.score.n."+TYPE+"."+output_file+".tsv"

elif TOOL == "SPLICEAI":
    TYPE = "noN"
    d_score_tsv_f = "../../src_tools_evaluation/spliceai_result/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
    a_score_tsv_f = "../../src_tools_evaluation/spliceai_result/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
    n_score_tsv_f = "../../src_tools_evaluation/spliceai_result/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"


FIGURE_ROOT = "Figures/"
def main():
    target_idx = 601
    # Read data from a TSV file

    donors = []
    with open(a_score_tsv_f, "r") as fr:
        lines = fr.read().splitlines()
        for line in lines:
            eles = line.split(" ")
            print("len(eles): ", len(eles))
            donors.append(float(eles[len(eles)-200]))
            # donors.append(float(eles[target_idx-1]))
            # print(len(donors))
    donors = np.array(donors)
    donors = pd.Series(donors)
    print("donors: ", donors)


    os.makedirs(FIGURE_ROOT, exist_ok=True)

    # d_df = pd.read_csv(d_score_tsv_f, sep=' ', header=None, squeeze=True)
    # print(d_df)
    # donors = d_df.iloc[:, target_idx]  # 0-based index, so 199 corresponds to the 200th column
    # # column_600 = d_df.iloc[:, 600]  # 0-based index, so 199 corresponds to the 200th column
    # print("donors: ", donors)
    # print("column_600: ", column_600)

    # # Generate random data following a normal distribution
    # # mean = 0
    # # std = 1
    # # data = np.random.normal(mean, std, 1000)  # generates 1000 random values from a normal distribution

    # Create a distribution plot with density plot
    sns.histplot(donors, kde=True, bins=30, color='blue', alpha=0.7)
    plt.xlabel('X-axis Label')
    plt.ylabel('Y-axis Label')
    plt.title('Distribution Plot with Normal Distribution')
    plt.grid(True)
    plt.savefig(FIGURE_ROOT+TOOL+"_"+output_file+"_acceptor.png")
    plt.show()

if __name__ == "__main__":
    main()