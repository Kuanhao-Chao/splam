import os
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text

annotations = ["chess", "gencode", "refseq_ucsc"]

colors = ['blue', 'red', 'green', 'cyan', 'magenta', 'olive', 'black', 'gray', 'orange', 'purple']

os.makedirs("exon_intron", exist_ok=True)

for level in ["exons", "introns", "loci"]:
    for library in ["polyA", "ribozero"]:
        for annotation in annotations:
            before_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/BEFORE.tsv", delimiter="\t", index_col=0)
            after_df = pd.read_csv("../../results/"+library+"/assembly/" + annotation + "/AFTER.tsv", delimiter="\t", index_col=0)

            rowname = list (after_df.index)
            print("rowname: ", rowname)

            texts = []
            for sample_id in range(10):

                a_df = after_df.iloc[[sample_id]]
                a_missed = a_df["missed_"+level]
                a_novel = a_df["novel_" + level]

                b_df = before_df.iloc[[sample_id]] 
                b_missed = b_df["missed_"+level]
                b_novel = b_df["novel_"+level]



                # # Existing data points
                # existing_Missed = [0.8, 0.7, 0.6, 0.5, 0.4]
                # existing_Novel = [0.9, 0.75, 0.65, 0.55, 0.35]

                # # New data point
                # new_Missed = 0.6
                # new_Novel = 0.8

                # Plotting the existing data points
                plt.scatter(b_novel, b_missed, color=colors[sample_id])

                # Plotting the new data point
                plt.scatter(a_novel, a_missed, color=colors[sample_id], label=rowname[sample_id])

                # Adding labels and title
                plt.axis('equal')
                plt.xlabel('Novel')
                plt.ylabel('Missed')
                plt.title('Missed vs. Novel (' + annotation + ')')

                # Adding an arrow between two dots
                arrow_start = (b_novel, b_missed)
                arrow_end = (a_novel, a_missed)
                plt.annotate("", xy=arrow_end, xytext=arrow_start, arrowprops=dict(arrowstyle='->', color=colors[sample_id]))

            plt.legend()
            plt.tight_layout()
            plt.savefig("exon_intron/"+library+ "_" + annotation + "_" + level + ".png", dpi=300)
            plt.close()
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

