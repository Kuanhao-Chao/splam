# Step 1: extract introns in the annotation
splam extract mouse_chr19.gff -o tmp_out_generalization

# Step 2: score introns in the annotation
splam score -A GRCm39_assembly_report.txt -G mouse_chr19.fa -m ../model/splam_script.pt -o tmp_out_generalization tmp_out_generalization/junction.bed

#Step 3: output statistics of each transcript
splam clean -o tmp_out_generalization
