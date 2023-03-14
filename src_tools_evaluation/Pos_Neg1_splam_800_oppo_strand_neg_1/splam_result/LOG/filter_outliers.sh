awk '{if ($5 == 1) {if ($5-$7 > 0.8 || $5-$8 > 0.8) {print} } else if($5 == 0) {if ($7-$5 > 0.8 || $8-$5 > 0.8) {print} }  } ' noshuffle_junction_score.bed
