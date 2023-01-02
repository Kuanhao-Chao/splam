SAMPLE=$1

samtools view /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.cleaned.bam | grep -o "NH:i:[0-9]*" > /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.cleaned.nh.txt

samtools view /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.discard.bam | grep -o "NH:i:[0-9]*" > /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/results/800bp/$SAMPLE/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v21/BAM/$SAMPLE.discard.nh.txt
