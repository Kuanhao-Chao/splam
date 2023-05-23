for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo samtools sort ../results/polyA/$sample/cleaned_2stage.bam  -o ../results/polyA/$sample/cleaned_2stage.sort.bam &
    samtools sort ../results/polyA/$sample/cleaned_2stage.bam  -o ../results/polyA/$sample/cleaned_2stage.sort.bam &
done
