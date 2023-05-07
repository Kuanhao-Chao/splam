awk '{if ($10 > 80 && ($15 > 50 && $16 > 50)) {print}}' show_coords.txt > show_coords_paralogs.txt
