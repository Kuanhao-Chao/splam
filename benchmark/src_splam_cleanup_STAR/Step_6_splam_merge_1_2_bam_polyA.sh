for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo samtools merge ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/$sample.cleaned.bam ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/cleaned.sort.bam ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/cleaned_2stage.sort.bam &
    samtools merge ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/$sample.cleaned.bam ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/cleaned.sort.bam ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/cleaned_2stage.sort.bam &
done
