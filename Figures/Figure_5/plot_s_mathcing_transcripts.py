import os
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

annotations = ["chess", "gencode", "refseq_ucsc"]

colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'olive', 'black', 'gray', 'orange', 'purple']

os.makedirs("sensitivity_mathcing_transcript", exist_ok=True)
for level in ["transcript", "locus"]:
    for library in ["polyA", "ribozero"]:
        for annotation in annotations:
            before_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/BEFORE.tsv", delimiter="\t", index_col=0)
            after_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/AFTER.tsv", delimiter="\t", index_col=0)

            rowname = list (after_df.index)
            print("rowname: ", rowname)
            plt.figure(figsize=(6, 6))  # Adjust the width and height as desired

            for sample_id in range(10):

                a_df = after_df.iloc[[sample_id]]
                # a_sensitivity = a_df[level+"_s"]
                if level == "transcript":
                    a_sensitivity = a_df["matching_transcripts"]
                else:
                    a_sensitivity = a_df["matching_loci"]

                a_precision = a_df[level+"_p"]

                b_df = before_df.iloc[[sample_id]] 
                # b_sensitivity = b_df[level+"_s"]
                if level == "transcript":
                    b_sensitivity = b_df["matching_transcripts"]
                else:
                    b_sensitivity = b_df["matching_loci"]

                b_precision = b_df[level+"_p"]

                # # Existing data points
                # existing_precision = [0.8, 0.7, 0.6, 0.5, 0.4]
                # existing_recall = [0.9, 0.75, 0.65, 0.55, 0.35]

                # # New data point
                # new_precision = 0.6
                # new_recall = 0.8

                # Plotting the existing data points
                plt.scatter(b_precision, b_sensitivity, color=colors[sample_id])

                # Plotting the new data point
                plt.scatter(a_precision, a_sensitivity, color=colors[sample_id], label=rowname[sample_id])

                # Adding labels and title
                # plt.axis('equal')
                plt.xlabel('Transcript precision (%)', labelpad=10)
                plt.ylabel('Assembled matching transcripts', labelpad=10)
                # # Set the pad parameter to increase the space between labels and ticks
                # plt.tick_params(axis='x', which='both', pad=10)
                # plt.tick_params(axis='y', which='both', pad=10)

                # plt.title('Precision vs. Recall (' + annotation + ')')
                if library == "polyA":
                    plt.title('Poly-A capture', fontsize = 20)
                elif library == "ribozero":
                    plt.title('Ribosomal RNA depletion', fontsize = 20)


                # Adding an arrow between two dots
                arrow_start = (b_precision, b_sensitivity)
                arrow_end = (a_precision, a_sensitivity)
                plt.annotate("", xy=arrow_end, xytext=arrow_start, arrowprops=dict(arrowstyle='->', color=colors[sample_id], mutation_scale=20))


            #     # Calculating the precision and recall difference
            #     precision_diff = float(a_precision - b_precision)
            #     recall_diff = float(a_sensitivity - b_sensitivity)

            #     # Displaying the precision and recall difference on top of the arrow
            #     arrow_mid = ((arrow_start[0] + arrow_end[0]) / 2, (arrow_start[1] + arrow_end[1]) / 2)

            #     # texts = [
            #     x = float((arrow_start[0] + arrow_end[0]) / 2)
            #     y = float((arrow_start[1] + arrow_end[1]) / 2)

            #     text = plt.text(x, y, f"({precision_diff:.2f}, {recall_diff:.2f})", ha='center', va='center')
            #     print("text: ", text)
            #     texts.append(text)

            # # text111 = [plt.text(1, 1, 'Text%s' %i, ha='center', va='center') for i in range(20)]
            # # print(text111)

            # print(texts)
            # adjust_text(texts)
            #     # plt.annotate(f"({precision_diff:.2f}, {recall_diff:.2f})", xy=text_coords, ha='center')


                # Adding a legend
            # plt.legend(bbox_to_anchor=(1.04, 0.96), loc="upper left")
            plt.tight_layout()
            plt.savefig("sensitivity_mathcing_transcript/"+library+ "_" + annotation + "_" + level + ".png", dpi=300)
            plt.close()