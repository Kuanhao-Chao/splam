# polyA
for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    # Stringtie assembly for before
    echo samtools index ../Dataset/polyA/$sample/$sample.bam &
    samtools index ../Dataset/polyA/$sample/$sample.bam &
done

# ribozero
for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287; do
    # Stringtie assembly for before
    echo samtools index ../Dataset/ribozero/$sample/$sample.bam &
    samtools index ../Dataset/ribozero/$sample/$sample.bam &
done