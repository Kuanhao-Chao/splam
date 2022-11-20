threshold=10
awk -v bound="$threshold" '{if($5 > bound && $6 != "?") {print}}' junctions_bed/junctions.bed > junctions_bed/junctions_$threshold.bed
