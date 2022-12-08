threshold=100
seq_len=1000
awk -v bound="$threshold" '{if($5 > bound && $6 != "?") {print}}' ../junctions_bed/junctions.bed > ../BAM_junctions/${seq_len}bp/junctions_$threshold.bed
