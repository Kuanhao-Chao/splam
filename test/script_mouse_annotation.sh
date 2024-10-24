# Step 1: extract introns in the annotation
splam extract mouse_chr19_subset.gff -o tmp_out_generalization

# Step 2: score introns in the annotation
splam score -G mouse_chr19.fa -m ../model/splam_static.pt -o tmp_out_generalization tmp_out_generalization/junction.bed

#Step 3: output statistics of each transcript
splam clean -o tmp_out_generalization -t 0.8
