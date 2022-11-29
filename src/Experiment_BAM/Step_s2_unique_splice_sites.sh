# awk '{print $1, $2, $3, "JUNC", "1", $6}' file.introns.bed | sort -k1,1 -k2,2n -k3,3n | uniq -u > introns.sort.bed

awk '{print $1, $2, $3, "JUNC", "1", $6}' SRR1352129_chr22.spliced.bed | sort -k1,1 -k2,2n -k3,3n | uniq -u > introns.sort.bed