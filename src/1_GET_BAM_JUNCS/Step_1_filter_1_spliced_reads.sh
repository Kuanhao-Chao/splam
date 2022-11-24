threshold=1
awk -v bound="$threshold" '{if($5 == bound && $6 != "?") {print}}' ../junctions_bed/junctions.bed > ../BAM_junctions/junctions_$threshold.bed
