import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns
import numpy as np
from collections import Counter

def main(argv):

    THRESHOLD = 0.2
    num_ls_discard = []
    num_ls_cleaned = []

    for sample in argv:
        print("sample: ", sample)
        
        j_score_bed_f = "../src_SPLAM/build/"+sample+"/junction_score.bed"
        # print("j_score_bed_f: ", j_score_bed_f)

        # discard_file="/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/"+sample+"/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/"+sample+".discard.juncs.bed"

        # cleaned_file="/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/"+sample+"/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/"+sample+".cleaned.juncs.bed"


        with open(j_score_bed_f, "r") as f:
            lines = f.read().splitlines()
            for line in lines:
                eles = line.split("\t")
                # print(eles)
                if len(eles) > 4:
                    num = int(eles[4])
                    score = float(eles[6])
                    # print(num)
                    # num = int(line[5:])
                    if score >= THRESHOLD:
                        num_ls_discard.append(num)
                    else :
                        num_ls_cleaned.append(num)
    # print("num_ls_discard: ", num_ls_discard)
    # print("num_ls_cleaned: ", num_ls_cleaned)
    num_ls_discard = np.array(num_ls_discard)
    num_ls_discard_counter = Counter(num_ls_discard)
    print("num_ls_discard_counter: ", num_ls_discard_counter)
    # num_ls_cleaned = np.array(num_ls_cleaned)
    # num_ls_cleaned_counter = Counter(num_ls_cleaned)

    keys_ls_discard = np.array(list(num_ls_discard_counter.keys()))
    vals_ls_discard = np.array(list(num_ls_discard_counter.values()))
    idices_discard = np.argsort(keys_ls_discard, axis=0)
    keys_ls_discard = keys_ls_discard[idices_discard]
    vals_ls_discard = vals_ls_discard[idices_discard]
    print("keys_ls discard: ", keys_ls_discard)
    print("vals_ls discard: ", vals_ls_discard)
    print("idices discard : ", idices_discard)

    # keys_ls_cleaned = np.array(list(num_ls_cleaned_counter.keys()))
    # vals_ls_cleaned = np.array(list(num_ls_cleaned_counter.values()))
    # idices_cleaned = np.argsort(keys_ls_cleaned, axis=0)
    # keys_ls_cleaned = keys_ls_cleaned[idices_cleaned]
    # vals_ls_cleaned = vals_ls_cleaned[idices_cleaned]
    # print("keys_ls cleaned: ", keys_ls_cleaned)
    # print("vals_ls cleaned: ", vals_ls_cleaned)
    # print("idices cleaned : ", idices_cleaned)

    # plt.bar(keys_ls_discard[:100]-0.2, vals_ls_discard[:100], width=0.2)
    # plt.bar(keys_ls_cleaned[:100]+0.2, vals_ls_cleaned[:100], width=0.2)

    # plt.bar(keys_ls_discard[:1000], vals_ls_discard[:1000])
    plt.bar(keys_ls_discard[:100], vals_ls_discard[:100])

    # plt.bar(keys_ls_cleaned[100:200]+0.2, vals_ls_cleaned[100:200], width=0.2)

    # plt.bar(num_ls_cleaned_counter.keys()[10:100], num_ls_cleaned_counter.values()[10:100], width=0.2)
    plt.xlabel("The alignment count for junctions")
    plt.ylabel("Ratio to 1-alignment junction")
    sample_name = "_".join(argv)
    plt.show()
    # plt.savefig(sample_name+".ratio.png")
    # plt.close()

    # # count_1_discard = len(num_ls_discard[num_ls_discard == 1])
    # # count_1_cleaned = len(num_ls_cleaned[num_ls_cleaned == 1])

    # # count_1_discard[count_1_discard > 100] = 100
    # count_1_discard = [x for x in num_ls_discard if x <= 100]
    # count_1_cleaned = [x for x in num_ls_cleaned if x <= 100]

    # ratio_ls_discard = []
    # # print("count_1: ", count_1)
    # for i in range(2,101):
    #     count_N_discard = len(num_ls_discard[num_ls_discard == i])
    #     ratio_ls_discard.append(count_N_discard/len(count_1_discard))
    #     print("count_"+str(i)+": ", count_N_discard)

    # ratio_ls_cleaned = []
    # # print("count_1: ", count_1)
    # for i in range(2,101):
    #     count_N_cleaned = len(num_ls_cleaned[num_ls_cleaned == i])
    #     ratio_ls_cleaned.append(count_N_cleaned/len(count_1_cleaned))
    #     print("count_"+str(i)+": ", count_N_cleaned)

    # plt.bar(np.arange(2,101)+0.2, ratio_ls_cleaned, width=0.2)
    # plt.bar(np.arange(2,101)-0.2, ratio_ls_discard, width=0.2)
    # plt.xlabel("The alignment count for junctions")
    # plt.ylabel("Ratio to 1-alignment junction")
    # sample_name = "_".join(argv)
    # plt.show()
    # plt.savefig(sample_name+".ratio.png")
    # plt.close()

    # ratio_ls = np.array(ratio_ls)
    # ratio_ls = np.insert(ratio_ls, 0, 1)
    # plt.bar([*range(1,101)], ratio_ls)
    # plt.xlabel("The alignment count for junctions")
    # plt.ylabel("Ratio to 1-alignment junction")
    # sample_name = "_".join(argv)
    # plt.savefig(sample_name+".ratio.with1.png")
    # plt.close()


    # sns.kdeplot(num_ls, label='Dicarded junctions')
    # # , clip=(0, 500))
    # # sns.kdeplot(num_ls_cleaned, label='Cleaned junctions')
    # # , clip = (0, 500))

    # sample_name = "_".join(argv)
    # plt.legend()
    # plt.savefig(sample_name+".all.png")
    # plt.close()
    # # plt.show()

    # num_ls[num_ls > 100] = 100
    # sns.kdeplot(num_ls, label='Dicarded junctions')
    # # , clip=(0, 500))
    # # sns.kdeplot(num_ls_cleaned, label='Cleaned junctions')
    # # , clip = (0, 500))

    # sample_name = "_".join(argv)
    # plt.legend()
    # plt.savefig(sample_name+".100.merged.all.png")
    # plt.close()
    # # plt.show()



if __name__ == "__main__":
    main(sys.argv[1:])

# samtools view /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/$SAMPLE.cleaned.bam | grep -o "NH:i:[0-9]*" > 

# samtools view /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/$SAMPLE.discard.bam | grep -o "NH:i:[0-9]*" > /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/$SAMPLE.discard.nh.txt
