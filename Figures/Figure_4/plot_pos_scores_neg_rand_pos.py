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
    output_files = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"] 
    FIGURE_ROOT = "Figures/"
    for SPLAM_VERSION in ["SPLAM_v11", "SPLAM_v12"]:
        os.makedirs(FIGURE_ROOT+SPLAM_VERSION+"/neg/", exist_ok=True)
        for output_file in output_files:
            fig = plt.figure(figsize=(8, 4))
            for TARGET in TARGETS:
                HANDELS = []
                for INDEX in range(len(TOOLS)):
                    TOOL = TOOLS[INDEX]
                    if TOOL == "SPLAM":
                        TYPE = "noshuffle"
                        d_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                        a_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                        n_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.n."+TYPE+"."+output_file+".tsv"
                    elif TOOL == "SPLICEAI":
                        TYPE = "noN"
                        d_score_tsv_f = "../../src_tools_evaluation/spliceai_result/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                        a_score_tsv_f = "../../src_tools_evaluation/spliceai_result/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                        n_score_tsv_f = "../../src_tools_evaluation/spliceai_result/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"


                    if TARGET == "Donor":
                        if TOOL == "SPLAM":
                            target_idx = 201
                        elif TOOL == "SPLICEAI":
                            target_idx = 200                        
                        donors = []


                        with open(d_score_tsv_f, "r") as fr:
                            lines = fr.read().splitlines()
                            counter = 0 
                            for line in lines:
                                counter += 1
                                if counter > 3000: 
                                    break
                                eles = line.split(" ")
                                print("len(eles): ", len(eles))

                                del eles[target_idx]
                                if TOOL == "SPLICEAI":
                                    eles = eles[0:399] + eles[len(eles)-400:len(eles)]
                                print(len(eles))

                                donors = donors + eles
                                # donors.append(float(eles[target_idx-1]))
                                # print(len(donors))
                        donors = np.array(donors)
                        donors = donors.astype(float)
                        print("donors: ", donors)
                        # Create a distribution plot with density plot
                        # sns.histplot(donors, kde=True, bins=50, color=COLORS[INDEX], alpha=0.7, stat="probability", label=TOOL)
                        # plt.hist(donors, density=True, histtype='stepfilled', alpha=0.8, label=TOOL)
                        
                        sns.kdeplot(donors, shade=True, clip = (0.0, 1.0), alpha=0.5, label=TOOL)

                    elif TARGET == "Acceptor":
                        target_idx = 200
                        acceptors = []

                        with open(a_score_tsv_f, "r") as fr:
                            lines = fr.read().splitlines()
                            counter = 0 
                            for line in lines:
                                counter += 1
                                if counter > 3000: 
                                    break
                                eles = line.split(" ")
                                print("len(eles): ", len(eles))
                                del eles[len(eles)-target_idx]

                                if TOOL == "SPLICEAI":
                                    eles = eles[0:400] + eles[len(eles)-399:len(eles)]
                                print(len(eles))
                                acceptors = acceptors + eles
                        acceptors = np.array(acceptors)
                        acceptors = acceptors.astype(float)
                        print("acceptors: ", acceptors)
                        # Create a distribution plot with density plot
                        # sns.histplot(acceptors, kde=True, bins=50, color=COLORS[INDEX], alpha=0.7, stat="probability", label=TOOL)
                        # plt.hist(acceptors, density=True, histtype='stepfilled', alpha=0.8, label=TOOL)
                        sns.kdeplot(acceptors, shade=True, clip = (0.0, 1.0), alpha=0.5, label=TOOL)
                    # HANDELS.append(plt_res)
                plt.legend(loc="upper center")
                plt.xlabel('Scores')
                plt.ylabel('Density')
                plt.title('Distribution Plot of '+TARGET+' Scores ('+output_file+')')
                plt.tight_layout()
                plt.grid(True)
                # Add a legend
                plt.savefig(FIGURE_ROOT+SPLAM_VERSION+"/neg/"+output_file+"_"+TARGET+".png", dpi=300)
                plt.clf()
                print(">>> Finish plotting 1 figure!")
                # plt.show()

if __name__ == "__main__":
    main()