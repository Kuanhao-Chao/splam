samtools view -h ../Dataset/SRR1352706_chr22.bam | awk '$0 ~ /^@/ || $6 ~ /N/' | samtools view -b > filtered.bam
