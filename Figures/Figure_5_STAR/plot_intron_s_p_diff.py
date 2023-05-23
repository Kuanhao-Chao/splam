import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

annotations = ["chess", "gencode", "refseq_ucsc"]

colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'olive', 'black', 'gray', 'orange', 'purple']

os.makedirs("intron_s_p_diff", exist_ok=True)
ind = np.array(list(range(10)))
width = 0.3
for library in ["polyA", "ribozero"]:
    for annotation in annotations:
        before_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/BEFORE.tsv", delimiter="\t", index_col=0)
        after_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/AFTER.tsv", delimiter="\t", index_col=0)

        rowname = list (after_df.index)
        print("rowname: ", rowname)

        difference_precision_ls = []
        difference_sensitivity_ls = []
        # plt.figure(figsize=(10, 3.5))  # Adjust the width and height as desired
        for sample_id in range(10):
            sample = rowname[sample_id]
            fr = open("../../results/"+library+"/"+sample+"/intron_matcher/BEFORE/res.txt", "r")
            b_line = fr.readline().splitlines()[0]
            b_eles = b_line.split("\t")
            print("b_eles: ", b_eles)
            b_precision = float(b_eles[4])
            b_sensitivity = float(b_eles[5])
            # for ele in b_eles:
            #     print(ele)

            fr = open("../../results/"+library+"/"+sample+"/intron_matcher/AFTER/res.txt", "r")
            a_line = fr.readline().splitlines()[0]
            a_eles = a_line.split("\t")
            print("a_eles: ", a_eles)
            a_precision = float(a_eles[4])
            a_sensitivity = float(a_eles[5])

            difference_precision_ls.append(a_precision - b_precision)
            difference_sensitivity_ls.append(a_sensitivity - b_sensitivity)


        difference_precision_ls = 100*np.array(difference_precision_ls)
        difference_sensitivity_ls = 100*np.array(difference_sensitivity_ls)

        colors_p = ['blue' if value > 0 else 'red' for value in difference_precision_ls]
        colors_s = ['blue' if value > 0 else 'red' for value in difference_sensitivity_ls]

        # Plotting the existing data points
        plt.figure(figsize=(10,6))
        plt.bar(ind, difference_precision_ls, width, color = colors_p, label="Precision difference (after - before)")
        plt.bar(ind + 0.3, difference_sensitivity_ls, width, color = colors_s, label="Sensitivity difference (after - before)")
        plt.xticks(ind + width / 2, rowname)
        plt.legend()
        # plt.xlabel('Samples')
        plt.ylabel("Difference(%)")
        if library == "polyA":
            plt.title("Poly-A capture", fontsize = 20)
        elif library == "ribozero":
            plt.title("Ribosomal RNA depletion", fontsize = 20)

        plt.tight_layout()
        plt.savefig("intron_s_p_diff/" + library + "_" + annotation + ".png", dpi=300)
