up_threshold=5
low_threshold=2
SEQ_LEN=800
mkdir BAM_junctions

awk -v up_bound="$up_threshold" -v low_bound="$low_threshold" '{if($5 <= up_bound && $5 >= low_bound && $6 != "?") {print}}' ../ALL_JUNCS/junctions.bed > ./BAM_junctions/junctions_${up_threshold}_${low_threshold}.bed
