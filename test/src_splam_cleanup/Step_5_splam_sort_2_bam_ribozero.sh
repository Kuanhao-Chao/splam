for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287; do
    echo samtools sort ../results/ribozero/$sample/cleaned.bam -o ../results/ribozero/$sample/cleaned.sort.bam &
    samtools sort ../results/ribozero/$sample/cleaned.bam -o ../results/ribozero/$sample/cleaned.sort.bam &
done
