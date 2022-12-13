# samtools view -h ../../Dataset/SRR1352129_chr22.bam  |  awk '($6 ~ /N/)' | samtools view -S -b > SRR1352129_chr22.spliced.bam

samtools view -h ../../Dataset/SRR1352129_chr22.bam | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b > SRR1352129_chr22.spliced.bam

# samtools view -S -b SRR1352129_chr22.spliced.sam > SRR1352129_chr22.spliced.bam