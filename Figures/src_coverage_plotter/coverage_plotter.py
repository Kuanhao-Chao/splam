import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def main():
    # with open("/Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/SRR1352129/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/SRR1352129.discard.pup", "r") as f:
    #     lines = f.read().splitlines()
    df_discard = pd.read_csv("/Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/SRR1352129/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/SRR1352129.discard.pup", delimiter='\t', header=None)

    df_clean = pd.read_csv("/Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/SRR1352129/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/BAM/SRR1352129.cleaned.pup", delimiter='\t', header=None)

    df_discard_1 = df_discard[df_discard[3] < 200]
    df_discard_200 = df_discard[df_discard[3] >= 200]
    
    df_clean_1 = df_clean[df_clean[3] < 200]
    df_clean_200 = df_clean[df_clean[3] >= 200]

    # print()
    
    # df2[df[3] < ]

    sns.kdeplot(df_discard_1[3])
    sns.kdeplot(df_clean_1[3])
    plt.savefig("./test_1.png", dpi=300)
    plt.close()

    sns.kdeplot(df_discard_200[3])
    sns.kdeplot(df_discard_200[3])
    plt.savefig("./test_2.png", dpi=300)
    plt.close()
    
    # sns.kdeplot(df2[3])
    # plt.plot(df[1], df[3])
    # plt.plot(df2[1], df2[3])

    # print(df)


if __name__ == "__main__":
    main()