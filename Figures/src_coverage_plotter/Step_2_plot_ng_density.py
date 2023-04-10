import matplotlib.pyplot as plt
import sys
import seaborn as sns

def main(argv):

    sample = argv[0]
    file="/Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/"+sample+"/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/"+sample+".discard.nh.txt"

    cleaned_file="/Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/"+sample+"/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/"+sample+".cleaned.nh.txt"

    num_ls = []
    with open(file, "r") as f:
        lines = f.read().splitlines()

        for line in lines:
            num = int(line[5:])
            num_ls.append(num)

    num_ls_cleaned = []
    with open(cleaned_file, "r") as f:
        lines = f.read().splitlines()

        for line in lines:
            num = int(line[5:])
            num_ls_cleaned.append(num)
    
    
    sns.kdeplot(num_ls)
    sns.kdeplot(num_ls_cleaned)

    plt.show()

if __name__ == "__main__":
    main(sys.argv[1:])

# samtools view /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.cleaned.bam | grep -o "NH:i:[0-9]*" > 

# samtools view /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.discard.bam | grep -o "NH:i:[0-9]*" > /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.discard.nh.txt
