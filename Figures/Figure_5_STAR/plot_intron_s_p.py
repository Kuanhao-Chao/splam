import os
import pandas as pd
import matplotlib.pyplot as plt

annotations = ["refseq_ucsc"]#["chess", "gencode", "refseq_ucsc"]

colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'olive', 'black', 'gray', 'orange', 'purple']

os.makedirs("intron_sensitivity_precision", exist_ok=True)
for library in ["polyA_STAR", "ribozero_STAR"]:
    for annotation in annotations:
        before_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/BEFORE.tsv", delimiter="\t", index_col=0)
        after_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/AFTER.tsv", delimiter="\t", index_col=0)

        rowname = list (after_df.index)
        print("rowname: ", rowname)
        plt.figure(figsize=(6, 2.57))  # Adjust the width and height as desired

        for sample_id in range(10):
            sample = rowname[sample_id]
            fr = open("../../results/"+library+"/"+sample+".bamAligned.sortedByCoord.out/intron_matcher/BEFORE/res.txt", "r")
            b_line = fr.readline().splitlines()[0]
            b_eles = b_line.split("\t")
            print("b_eles: ", b_eles)
            b_precision = 100*float(b_eles[4])
            b_sensitivity = 100*float(b_eles[5])
            # for ele in b_eles:
            #     print(ele)

            fr = open("../../results/"+library+"/"+sample+".bamAligned.sortedByCoord.out/intron_matcher/AFTER/res.txt", "r")
            a_line = fr.readline().splitlines()[0]
            a_eles = a_line.split("\t")
            print("a_eles: ", a_eles)
            a_precision = 100*float(a_eles[4])
            a_sensitivity = 100*float(a_eles[5])
        
            

            # Plotting the existing data points
            plt.scatter(b_precision, b_sensitivity, color=colors[sample_id])

            # Plotting the new data point
            plt.scatter(a_precision, a_sensitivity, color=colors[sample_id], label=rowname[sample_id])

            # Adding labels and title
            # plt.rcParams["figure.autolayout"] = True
            plt.axis('equal')
            # plt.ylim(50, 60)
            plt.xlabel('Intron precision (%)', labelpad=10)
            plt.ylabel('Intron recall (%)', labelpad=10)
            # plt.title('Precision vs. Recall (' + annotation + ')')
            if library == "polyA":
                plt.title("Poly-A capture", fontsize = 18)
            elif library == "ribozero":
                plt.title("Ribosomal RNA depletion", fontsize = 18)

            # Adding an arrow between two dots
            arrow_start = (b_precision, b_sensitivity)
            arrow_end = (a_precision, a_sensitivity)
            # plt.annotate("", xy=arrow_end, xytext=arrow_start, arrowprops=dict(arrowstyle='->', color=colors[sample_id]))
            plt.annotate("", xy=arrow_end, xytext=arrow_start, arrowprops=dict(arrowstyle='->', color=colors[sample_id], mutation_scale=20))

        # plt.legend()
        # plt.legend(bbox_to_anchor=(1.04, 0.96), loc="upper left")
        # plt.legend(bbox_to_anchor=(0.5, 1.25), loc="upper center", ncol=5)
        plt.tight_layout()
        plt.savefig("intron_sensitivity_precision/"+library+ "_" + annotation + ".png", dpi=300)
        plt.close()
            # Displaying the plot
        # plt.show()
