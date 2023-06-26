samtools view -h ../../Dataset/SRR1352415.bam | awk '$0 ~ /^@/ || $6 ~ /N/'  \
   | samtools view -bS - \
   | bamToBed -bed12 > file.12.bed

bed12ToBed6 -i file.12.bed > file.6.bed
subtractBed -a file.12.bed -b file.6.bed -s | cut -f 1-6 > file.introns.bed