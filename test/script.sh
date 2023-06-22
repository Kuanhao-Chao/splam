# Step 1
splam extract -P SRR1352129_chr9_sub.bam

# Step 2
splam score -G chr9.fa -m ../model/splam_script.pt -o tmp_splam_out tmp_splam_out/junction.bed

#Step 3
splam clean -o tmp_splam_out