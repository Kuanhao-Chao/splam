threshold=1
SEQ_LEN=600
awk -v bound="$threshold" '{if($5 == bound && $6 != "?") {print}}' ../junctions_bed/junctions.bed > ../BAM_junctions/${SEQ_LEN}bp/junctions_$threshold.bed
