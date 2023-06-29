# polyA
for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo mkdir ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/
    mkdir ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/
    
    echo stringtie -o ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/$sample.gtf ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/$sample.cleaned.bam &
    stringtie -o ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/$sample.gtf ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/$sample.cleaned.bam &
done

# ribozero
for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287; do
    echo mkdir ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/
    mkdir ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/

    echo stringtie -o ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/$sample.gtf ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/$sample.cleaned.bam &
    stringtie -o ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/$sample.gtf ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/$sample.cleaned.bam &
done