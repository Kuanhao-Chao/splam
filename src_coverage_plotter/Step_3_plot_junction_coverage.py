import matplotlib.pyplot as plt
import pandas as pd
import sys
import seaborn as sns
import numpy as np

def main(argv):

    num_ls = []
    num_ls_cleaned = []

    for sample in argv:
        print("sample: ", sample)
        discard_file="/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/"+sample+"/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/"+sample+".discard.juncs.bed"

        cleaned_file="/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/"+sample+"/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/"+sample+".cleaned.juncs.bed"

        with open(discard_file, "r") as f:
            lines = f.read().splitlines()
            for line in lines:
                eles = line.split("\t")
                if len(eles) > 4:
                    num = int(eles[4])
                    # print(num)
                    # num = int(line[5:])
                    num_ls.append(num)

        with open(cleaned_file, "r") as f:
            lines = f.read().splitlines()
            for line in lines:
                eles = line.split("\t")
                if len(eles) > 4:
                    num = int(eles[4])
                    # print(num)
                    # num = int(line[5:])
                    num_ls_cleaned.append(num)
    
    num_ls = np.array(num_ls)
    num_ls_cleaned = np.array(num_ls_cleaned)

    num_ls = num_ls[num_ls <= 100]
    num_ls_cleaned = num_ls_cleaned[num_ls_cleaned <= 100 ]

    # num_ls[num_ls > 100] = 0
    # num_ls_cleaned[num_ls_cleaned > 100 ] = 0

    # # num_ls_selected = num_ls[num_ls > 10]


    # # num_ls_1 = num_ls[num_ls <= 20]
    # # num_ls_200 = num_ls[num_ls > 200]
    # # # num_ls_selected = num_ls[num_ls > 200 & num_ls < 1000]
    # # num_ls_selected = num_ls[np.where((num_ls > 200) & (num_ls < 1000))]

    # # num_ls_cleaned_1 = num_ls_cleaned[num_ls_cleaned <= 20]
    # # num_ls_cleaned_200 = num_ls_cleaned[num_ls_cleaned > 200]
    # # # num_ls_cleaned_selected = num_ls_cleaned[num_ls_cleaned > 200 & num_ls_cleaned < 1000]
    # # num_ls_cleaned_selected = num_ls_cleaned[np.where((num_ls_cleaned > 200) & (num_ls_cleaned < 1000))]


    # # print(num_ls_cleaned_1)

    # # sns.distplot(num_ls, hist=True, kde=False)
    # # sns.distplot(num_ls_cleaned, hist=True, kde=False)

    #         #  bins=int(180/5), color = 'blue',
    #         #  hist_kws={'edgecolor':'black'})
    
    # sns.distplot(num_ls, hist=True, kde_kws={'clip': (0.0, 100)})
    # sns.distplot(num_ls_cleaned, hist=True, kde_kws={'clip': (0.0, 100)})

    sns.kdeplot(num_ls, label='Dicarded junctions')
    # , clip=(0, 500))
    sns.kdeplot(num_ls_cleaned, label='Cleaned junctions')
    # , clip = (0, 500))

    sample_name = "_".join(argv)
    plt.legend()
    plt.savefig(sample_name+".png")
    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])

# samtools view /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.cleaned.bam | grep -o "NH:i:[0-9]*" > 

# samtools view /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.discard.bam | grep -o "NH:i:[0-9]*" > /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.discard.nh.txt
