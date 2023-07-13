import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import random
import os

def main():
    # for TOOL in ["SPLAM", "SPLICEAI"]:
    COLORS = ["#2ca02c", "#ff7f0e"]
    TOOLS = ["splam", "SpliceAI"]#, "SpliceAI-10k-Ns"]
    TARGETS = ["Donor", "Acceptor"]
    output_files = ["pos_MANE", "pos_ALTS", "neg_1", "neg_random"] 
    # output_files = ["pos_MANE"]#, "neg_1", "neg_random"] 
    FIGURE_ROOT = "IMG/"
    SPLAM_VERSION = "SPLAM_v11"


    os.makedirs(FIGURE_ROOT+"AVERAGE/", exist_ok=True)

    for output_file in output_files:

        splam_donors_sum = np.zeros(10000)
        spliceai_donors_sum = np.zeros(10000)
        splam_acceptors_sum = np.zeros(10000)
        spliceai_acceptors_sum = np.zeros(10000)
        # for SPLICEAI_VERSION in ["1"]:
        for SPLICEAI_VERSION in ["1", "2", "3", "4", "5"]:
            fig = plt.figure(figsize=(6, 3))
            for TARGET in TARGETS:
                HANDELS = []
                for INDEX in range(len(TOOLS)):
                    TOOL = TOOLS[INDEX]
                    if TOOL == "splam":
                        TYPE = "noshuffle"
                        d_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                        a_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                        n_score_tsv_f = "../../src_tools_evaluation/splam_result/"+SPLAM_VERSION+"/"+output_file+"/splam_all_seq.score.n."+TYPE+"."+output_file+".tsv"


                    elif TOOL == "SpliceAI":
                        TYPE = "noN"
                        d_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                        a_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                        n_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"
                    
                    elif TOOL == "SpliceAI-10k-Ns":
                        TYPE = "N"
                        d_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.d."+TYPE+"."+output_file+".tsv"
                        a_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.a."+TYPE+"."+output_file+".tsv"
                        n_score_tsv_f = "../../src_tools_evaluation/spliceai_result_"+SPLICEAI_VERSION+"/"+output_file+"/spliceai_all_seq.score.n."+TYPE+"."+output_file+".tsv"

                    print("d_score_tsv_f: ", d_score_tsv_f)
                    print("a_score_tsv_f: ", a_score_tsv_f)
                    print("n_score_tsv_f: ", n_score_tsv_f)

                    if TARGET == "Donor":
                        if TOOL == "splam":
                            target_idx = 201
                        elif TOOL == "SpliceAI" or TOOL == "SpliceAI-10k-Ns":
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

                        if TOOL == "splam":
                            splam_donors_sum += donors
                        elif TOOL == "SpliceAI" or TOOL == "SpliceAI-10k-Ns":
                            spliceai_donors_sum += donors

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

                        if TOOL == "splam":
                            splam_acceptors_sum += acceptors
                        elif TOOL == "SpliceAI" or TOOL == "SpliceAI-10k-Ns":
                            spliceai_acceptors_sum += acceptors


        # Create a distribution plot with density plot
        splam_donors_sum /= 5
        spliceai_donors_sum /= 5
        splam_acceptors_sum /= 5
        spliceai_acceptors_sum /= 5


        for target_tool in ["splam", "SpliceAI"]:
            # if target_tool == "splam":
            #     d_score_ls = splam_donors_sum
            #     a_score_ls = splam_acceptors_sum

            # elif target_tool == "SpliceAI":
            #     d_score_ls = spliceai_donors_sum
            #     a_score_ls = spliceai_acceptors_sum

            color_splam = "#2ca02c" 
            color_spliceai = "#ff7f0e"


        # plot donor
        TARGET = "Donor"
        sns.kdeplot(splam_donors_sum, shade=True, clip = (0.0, 1.0), alpha=0.5, label="Splam", color=color_splam)
        sns.kdeplot(spliceai_donors_sum, shade=True, clip = (0.0, 1.0), alpha=0.5, label="SpliceAI", color=color_spliceai)
        plt.legend(loc="upper center")
        plt.xlabel('Scores')
        plt.ylabel('Density')
        plt.title('Distribution Plot of '+TARGET+' Scores ('+output_file+')')
        plt.tight_layout()
        plt.grid(True)
        # Add a legend
        plt.savefig(FIGURE_ROOT+"AVERAGE/"+output_file+"_"+TARGET+".png", dpi=300)
        plt.clf()
        print(">>> Finish plotting 1 figure!")
        # plt.show()

        # plot acceptor
        TARGET = "Acceptor"
        sns.kdeplot(splam_acceptors_sum, shade=True, clip = (0.0, 1.0), alpha=0.5, label="Splam", color=color_splam)
        sns.kdeplot(spliceai_acceptors_sum, shade=True, clip = (0.0, 1.0), alpha=0.5, label="SpliceAI", color=color_spliceai)
        plt.legend(loc="upper center")
        plt.xlabel('Scores')
        plt.ylabel('Density')
        plt.title('Distribution Plot of '+TARGET+' Scores ('+output_file+')')
        plt.tight_layout()
        plt.grid(True)
        # Add a legend
        plt.savefig(FIGURE_ROOT+"AVERAGE/"+output_file+"_"+TARGET+".png", dpi=300)
        plt.clf()
        print(">>> Finish plotting 1 figure!")
        # plt.show()



if __name__ == "__main__":
    main()