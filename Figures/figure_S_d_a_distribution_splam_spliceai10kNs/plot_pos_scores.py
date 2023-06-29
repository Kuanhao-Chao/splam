import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import random
import os

def main():
    # for TOOL in ["SPLAM", "SPLICEAI"]:
    COLORS = ["green", "blue"]
    TOOLS = ["SPLAM", "SpliceAI-10k-Ns"]#, "SpliceAI-10k-Ns"]
    TARGETS = ["Donor", "Acceptor"]
    output_files = ["pos_MANE", "pos_ALTS"]#, "neg_1", "neg_random"] 
    FIGURE_ROOT = "Figures/"
    for SPLICEAI_VERSION in ["1", "2", "3", "4", "5", "AVERAGE"]:
        for SPLAM_VERSION in ["SPLAM_v11"]:#, "SPLAM_v12"]:
            os.makedirs(FIGURE_ROOT+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/", exist_ok=True)
            for output_file in output_files:
                fig = plt.figure(figsize=(6, 4))
                for TARGET in TARGETS:
                    HANDELS = []
                    for INDEX in range(len(TOOLS)):
                        TOOL = TOOLS[INDEX]
                        if TOOL == "SPLAM":
                            TYPE = "noshuffle"
                            d_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                            a_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                            n_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.n."+TYPE+"."+output_file+".tsv"


                        elif TOOL == "SpliceAI-10k":
                            TYPE = "noN"
                            d_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                            a_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                            n_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"
                        
                        elif TOOL == "SpliceAI-10k-Ns":
                            TYPE = "N"
                            d_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                            a_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                            n_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"

                        color_plt = ""
                        if TOOL == "SPLAM":
                            color_plt = "#2ca02c" 
                        elif TOOL == "SpliceAI-10k" or TOOL == "SpliceAI-10k-Ns":
                            color_plt = "#1f77b4"

                        if TARGET == "Donor":
                            if TOOL == "SPLAM":
                                target_idx = 201
                            elif TOOL == "SpliceAI-10k" or TOOL == "SpliceAI-10k-Ns":
                                target_idx = 200
                            donors = []
                            with open(d_score_tsv_f, "r") as fr:
                                lines = fr.read().splitlines()
                                for line in lines:
                                    eles = line.split(" ")
                                    print("len(eles): ", len(eles))
                                    donors.append(float(eles[target_idx-1]))
                                    # print(len(donors))
                            donors = np.array(donors)
                            print("donors: ", donors)
                            # Create a distribution plot with density plot
                            # sns.histplot(donors, kde=True, bins=50, color=COLORS[INDEX], alpha=0.7, stat="probability", label=TOOL)
                            # plt.hist(donors, density=True, histtype='stepfilled', alpha=0.8, label=TOOL)
                            
                            sns.kdeplot(donors, shade=True, clip = (0.0, 1.0), alpha=0.5, label=TOOL, color=color_plt)

                        elif TARGET == "Acceptor":
                            target_idx = 200
                            acceptors = []
                            with open(a_score_tsv_f, "r") as fr:
                                lines = fr.read().splitlines()
                                for line in lines:
                                    eles = line.split(" ")
                                    acceptors.append(float(eles[len(eles)-target_idx]))
                                    # acceptors.append(float(eles[target_idx-1]))
                                    # print(len(acceptors))
                            acceptors = np.array(acceptors)
                            print("acceptors: ", acceptors)
                            # Create a distribution plot with density plot
                            # sns.histplot(acceptors, kde=True, bins=50, color=COLORS[INDEX], alpha=0.7, stat="probability", label=TOOL)
                            # plt.hist(acceptors, density=True, histtype='stepfilled', alpha=0.8, label=TOOL)
                            sns.kdeplot(acceptors, shade=True, clip = (0.0, 1.0), alpha=0.5, label=TOOL, color=color_plt)
                        # HANDELS.append(plt_res)
                    plt.legend(loc="upper center")
                    plt.xlabel('Scores')
                    plt.ylabel('Density')
                    plt.title('Distribution Plot of '+TARGET+' Scores ('+output_file+')')
                    plt.tight_layout()
                    plt.grid(True)
                    # Add a legend
                    plt.savefig(FIGURE_ROOT+SPLAM_VERSION+"/"+SPLICEAI_VERSION+"/" + output_file+"_"+TARGET+".png", dpi=300)
                    plt.clf()
                    print(">>> Finish plotting 1 figure!")
                    # plt.show()

if __name__ == "__main__":
    main()