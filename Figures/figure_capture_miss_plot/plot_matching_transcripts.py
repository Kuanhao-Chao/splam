import os
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

annotations = ["chess", "gencode", "refseq_ucsc"]

colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'olive', 'black', 'gray', 'orange', 'purple']

os.makedirs("matching_stats", exist_ok=True)
for level in ["transcripts", "loci"]:
    for library in ["polyA", "ribozero"]:
        for annotation in annotations:
            print("../../results/"+library+"/assembly/" + annotation + "/BEFORE.tsv")
            before_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/BEFORE.tsv", delimiter="\t", index_col=0)
            after_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/AFTER.tsv", delimiter="\t", index_col=0)

            rowname = list (after_df.index)
            print("rowname: ", rowname)

            texts = []
            a_matching_ls = []
            b_matching_ls = []
            difference_matching_ls = []
            for sample_id in range(10):

                a_df = after_df.iloc[[sample_id]]
                a_matching = float(a_df["matching_" + level])
                a_matching_ls.append(a_matching)

                b_df = before_df.iloc[[sample_id]] 
                b_matching = float(b_df["matching_" + level])
                b_matching_ls.append(b_matching)
                difference_matching_ls.append(100*(a_matching - b_matching) / b_matching)

            print("a_matching_ls: ", a_matching_ls)
            print("b_matching_ls: ", b_matching_ls)

            colors = ['blue' if value > 0 else 'red' for value in difference_matching_ls]
            # Plotting the existing data points
            plt.figure(figsize=(10,6))
            plt.bar(rowname, difference_matching_ls, color = colors)
            plt.xlabel('Samples')
            plt.ylabel('Number of matching ' + level + " (%)")
            plt.title('Precision vs. Recall (' + annotation + ')')

            plt.tight_layout()
            plt.savefig("matching_stats/" + level + "_" + library + "_" + annotation + ".png", dpi=300)

            # # Plotting the new data point
            # plt.scatter(a_precision, a_sensitivity, color=colors[sample_id], label=rowname[sample_id])

            #     # Adding labels and title
            #     plt.axis('equal')
            #     plt.xlabel('Precision')
            #     plt.ylabel('Recall')
            #     plt.title('Precision vs. Recall (' + annotation + ')')

            #     # Adding an arrow between two dots
            #     arrow_start = (b_precision, b_sensitivity)
            #     arrow_end = (a_precision, a_sensitivity)
            #     plt.annotate("", xy=arrow_end, xytext=arrow_start, arrowprops=dict(arrowstyle='->', color=colors[sample_id]))

            # #     # Calculating the precision and recall difference
            # #     precision_diff = float(a_precision - b_precision)
            # #     recall_diff = float(a_sensitivity - b_sensitivity)

            # #     # Displaying the precision and recall difference on top of the arrow
            # #     arrow_mid = ((arrow_start[0] + arrow_end[0]) / 2, (arrow_start[1] + arrow_end[1]) / 2)

            # #     # texts = [
            # #     x = float((arrow_start[0] + arrow_end[0]) / 2)
            # #     y = float((arrow_start[1] + arrow_end[1]) / 2)

            # #     text = plt.text(x, y, f"({precision_diff:.2f}, {recall_diff:.2f})", ha='center', va='center')
            # #     print("text: ", text)
            # #     texts.append(text)

            # # # text111 = [plt.text(1, 1, 'Text%s' %i, ha='center', va='center') for i in range(20)]
            # # # print(text111)

            # # print(texts)
            # # adjust_text(texts)
            # #     # plt.annotate(f"({precision_diff:.2f}, {recall_diff:.2f})", xy=text_coords, ha='center')


            #     # Adding a legend
            # plt.legend()
            # plt.tight_layout()
            # plt.savefig(library+ "_" + annotation + "_" + level + ".png", dpi=300)
            # plt.close()
                # Displaying the plot
            # plt.show()

                # print("a_df['transcript_s']: ", a_df["transcript_s"])
                # print("a_df['transcript_p']", a_df["transcript_p"])

    #             a_df = a_df.values.tolist()[0]
    #             a_df = [float(value) for value in a_df]
    #             print("a_df: ", a_df)


    #             # diff_df = after_df.iloc[[sample_id]] - before_df.iloc[[sample_id]] 
    #             # print(diff_df)
    #             # diff_df[""]
    #             # y = diff_df.values.tolist()[0]
    #             # y = [float(value) for value in y]
    #             # print("y: ", y)

    #             # Create a list of colors based on the values
    #             colors = ['blue' if value > 0 else 'red' for value in y]

    #             # Create the bar chart
    #             plt.bar(list(before_df.columns.values)[:12] , y[:12], color = colors)

    #             # Add labels and title
    #             plt.xlabel('X-axis')
    #             plt.ylabel('Y-axis')
    #             plt.title('Bar Chart with Blue and Red Bars')

    #             # Display the chart
    #             plt.show()
    #             # print("before_df: ", before_df.iloc[[4]])
    #             # print("after_df : ", after_df.iloc[[4]])

