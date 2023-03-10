threshold=1
SEQ_LEN=800
mkdir BAM_junctions

awk -v bound="$threshold" '{if($5 == bound && $6 != "?") {print}}' ../ALL_JUNCS/junctions.bed > ./BAM_junctions/junctions_$threshold.bed
