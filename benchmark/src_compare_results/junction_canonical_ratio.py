import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


for type in ["donor", "acceptor"]:
    targets = ["5_GET_REF_JUNCS_MANE", "6_GET_REF_JUNCS_REFSEQ_ALTS", "4_GET_BAM_NEG_OPP_STRAND_JUNCS"]
    X_axis = ["Positive_MANE", "Positive_ALT", "Negative_1"]
    barWidth = 0.25
    dataframes = []
    donor_top1 = []
    donor_top2 = []
    donor_top3 = []
    donor_top4 = []
    donor_top1_seq = []
    donor_top2_seq = []
    donor_top3_seq = []
    donor_top4_seq = []

    acceptor_top1 = []
    acceptor_top2 = []
    acceptor_top3 = []
    acceptor_top4 = []

    acceptor_top1_seq = []
    acceptor_top2_seq = []
    acceptor_top3_seq = []
    acceptor_top4_seq = []

    sum_grp = [0, 0, 0, 0]
    labels = ["Top1", "Top2", "Top3", "Top4"]

    idx = 0
    for target in targets:
        # print("target: ", target)
        # df = pd.read_csv(target+"/d_a_type.tsv", sep='\t', header=None)
        # dataframes.append(df)

    # column_names = ['Pair', 'Count']

    # # Create an empty DataFrame to store the combined data
    # combined_df = pd.DataFrame(columns=column_names)

    # # Load each file and add its data to the combined DataFrame
    # for target in targets:
    #     print("target: ", target)
    #     df = pd.read_csv(target+"/d_a_type.tsv", sep='\t', header=None)
    #     df.columns = column_names
    #     combined_df = pd.concat([combined_df, df])

    # print("combined_df: ", combined_df)
    # # Pivot the data to create a DataFrame with pairs as columns and counts as values
    # pivot_df = combined_df.pivot(index=None, columns='Pair', values='Count')

    # # Create a stacked bar chart of the data
    # pivot_df.plot(kind='bar', stacked=True)

    # # Add labels and a title to the chart
    # plt.xlabel('Sample')
    # plt.ylabel('Count')
    # plt.title('Counts of nucleotide pairs in three samples')

    # # Show the chart
    # plt.show()


        with open(target+"/d_a_type.tsv", "r") as fr:
            lines = fr.read().splitlines()
            donors = []
            donors_score = []

            if type == "donor":
                eles1 = lines[0].split("\t")
                donor_top1_seq.append(eles1[0])
                # donor_top1.append(int(eles1[1]))
                sum_grp[idx] += int(eles1[1])
                # /int(eles1[1]))


                eles2 = lines[1].split("\t")
                donor_top2_seq.append(eles2[0])
                # donor_top2.append(int(eles2[1]))
                sum_grp[idx] += int(eles2[1])
                # / int(eles2[1]))


                eles3 = lines[2].split("\t")
                donor_top3_seq.append(eles3[0])
                sum_grp[idx] += int(eles3[1])

                eles4 = lines[3].split("\t")
                donor_top4_seq.append(eles4[0])
                sum_grp[idx] += int(eles4[1])


                donor_top1.append(int(eles1[1]) / sum_grp[idx])
                donor_top2.append(int(eles2[1]) / sum_grp[idx])
                donor_top3.append(int(eles3[1]) / sum_grp[idx])
                donor_top4.append(int(eles4[1]) / sum_grp[idx])
                idx += 1
            elif type == "acceptor":
                eles1 = lines[5].split("\t")
                donor_top1_seq.append(eles1[0])
                # donor_top1.append(int(eles1[1]))
                sum_grp[idx] += int(eles1[1])
                # /int(eles1[1]))


                eles2 = lines[6].split("\t")
                donor_top2_seq.append(eles2[0])
                # donor_top2.append(int(eles2[1]))
                sum_grp[idx] += int(eles2[1])
                # / int(eles2[1]))


                eles3 = lines[7].split("\t")
                donor_top3_seq.append(eles3[0])
                sum_grp[idx] += int(eles3[1])

                eles4 = lines[8].split("\t")
                donor_top4_seq.append(eles4[0])
                sum_grp[idx] += int(eles4[1])


                donor_top1.append(int(eles1[1]) / sum_grp[idx])
                donor_top2.append(int(eles2[1]) / sum_grp[idx])
                donor_top3.append(int(eles3[1]) / sum_grp[idx])
                donor_top4.append(int(eles4[1]) / sum_grp[idx])
                idx += 1
    

    print("donor_top1: ", donor_top1)
    print("donor_top2: ", donor_top2)
    print("donor_top3: ", donor_top3)
    print("donor_top4: ", donor_top4)

    X = np.arange(3)
    fig, ax = plt.subplots(figsize=(10, 4))

    plt.ylim(0, 1.1)

    plt.bar(X-0.3, donor_top1, width = 0.2, label=labels[0])
    for i, total in enumerate(donor_top1):
        print(i, total)
        plt.text(i-0.3, donor_top1[i] + 0.01, donor_top1_seq[i] + "\n" + str(round(total*100, 2)) + "%", ha='center',  weight='bold')
        #ax.text(totals.index[i], total + y_offset, round(total), ha='center',  weight='bold')
    # plt.text(X-0.3, donor_top1, f'{donor_top1:.0f}', ha='center', va='bottom')
    # plt.text(0.5, 0.5, 'matplotlib', horizontalalignment='center')
        
    plt.bar(X-0.1, donor_top2, width = 0.2, label=labels[1])
    for i, total in enumerate(donor_top2):
        print(i, total)
        plt.text(i-0.1, donor_top2[i] + 0.01, donor_top2_seq[i] + "\n" + str(round(total*100, 2)) + "%", ha='center',  weight='bold')

    # , bottom=donor_top1)
    plt.bar(X+0.10, donor_top3, width = 0.2, label=labels[2])
    for i, total in enumerate(donor_top3):
        print(i, total)
        plt.text(i+0.1, donor_top3[i] + 0.01, donor_top3_seq[i] + "\n" + str(round(total*100, 2)) + "%", ha='center',  weight='bold')

    plt.bar(X+0.30, donor_top4, width = 0.2, label=labels[3])
    for i, total in enumerate(donor_top4):
        print(i, total)
        plt.text(i+0.3, donor_top4[i] + 0.01, donor_top4_seq[i] + "\n" + str(round(total*100, 2)) + "%", ha='center',  weight='bold')

    plt.xticks(X, X_axis)
    plt.legend()
    # xlabel(r'\textbf{X-AXIS}', fontsize=20)

    ax.xaxis.set_tick_params(labelsize=10, pad=10)
    ax.yaxis.set_tick_params(labelsize=10, pad=10)
    plt.ylabel("Ratio")
    plt.tight_layout()

    plt.savefig("junction_canonical_ratio_"+type+".png", dpi=300)
