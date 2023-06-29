awk '{if ($10 > 80 && ($16 > 50)) {print}}' show_coords.txt > show_coords_paralogs.txt
awk '{ print $19 }' show_coords_paralogs.txt > show_coords_paralogs.bed
