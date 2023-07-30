# Step 1: extract introns in the annotation
splam extract MANE.GRCh38.v1.1.subset.gff -o tmp_out_annotation

# Step 2: score introns in the annotation
splam score -G chr9_subset.fa -m ../model/splam_script.pt -o tmp_out_annotation tmp_out_annotation/junction.bed

#Step 3: output statistics of each transcript
splam clean -o tmp_out_annotation
