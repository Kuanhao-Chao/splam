import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import random
import os

def main():
    # for TOOL in ["SPLAM", "SPLICEAI"]:
    COLORS = ["green", "blue"]
    TOOLS = ["SPLAM", "SPLICEAI"]
    TARGETS = ["Donor", "Acceptor"]
    output_file = "pos_ALTS"
    fig, ax = plt.subplots()
    for TARGET in TARGETS:
        HANDELS = []
        for INDEX in range(len(TOOLS)):
            TOOL = TOOLS[INDEX]
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
            os.makedirs(FIGURE_ROOT, exist_ok=True)
            if TARGET == "Donor":
                if TOOL == "SPLAM":
                    target_idx = 201
                elif TOOL == "SPLICEAI":
                    target_idx = 200
                donors = []
                with open(d_score_tsv_f, "r") as fr:
                    lines = fr.read().splitlines()
                    for line in lines:
                        eles = line.split(" ")
                        print("len(eles): ", len(eles))

                        rand_idx = target_idx
                        while rand_idx != target_idx:
                            rand_idx = random.random()
                        donors.append(float(eles[target_idx-1]))
                        # print(len(donors))
                donors = np.array(donors)
                print("donors: ", donors)
                # Create a distribution plot with density plot
                # sns.histplot(donors, kde=True, bins=50, color=COLORS[INDEX], alpha=0.7, stat="probability", label=TOOL)
                # plt.hist(donors, density=True, histtype='stepfilled', alpha=0.8, label=TOOL)
                sns.kdeplot(donors, shade=True)


            elif TARGET == "Acceptor":
                target_idx = 200
                acceptors = []
                with open(a_score_tsv_f, "r") as fr:
                    lines = fr.read().splitlines()
                    for line in lines:
                        eles = line.split(" ")
                        print("len(eles): ", len(eles))
                        acceptors.append(float(eles[len(eles)-target_idx]))
                        # acceptors.append(float(eles[target_idx-1]))
                        # print(len(acceptors))
                acceptors = np.array(acceptors)
                print("acceptors: ", acceptors)
                # Create a distribution plot with density plot
                # sns.histplot(acceptors, kde=True, bins=50, color=COLORS[INDEX], alpha=0.7, stat="probability", label=TOOL)
                # plt.hist(acceptors, density=True, histtype='stepfilled', alpha=0.8, label=TOOL)
                sns.kdeplot(acceptors, shade=True)
            # HANDELS.append(plt_res)
        plt.legend()
        plt.xlabel('Scores')
        plt.ylabel('Probability')
        plt.title('Distribution Plot of '+TARGET+' Scores ('+output_file+')')
        plt.grid(True)
        # Add a legend
        plt.savefig(FIGURE_ROOT+output_file+"_"+TARGET+".png", dpi=300)
        plt.clf()
        print(">>> Finish plotting 1 figure!")
        # plt.show()

if __name__ == "__main__":
    main()