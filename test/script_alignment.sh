# Step 1: extract splice junctions in the alignment file
splam extract -P SRR1352129_chr9_sub.bam -o tmp_out_alignment --fr

# Step 2: score all the extracted splice junctions
splam score -G chr9_subset.fa -m ../model/splam_static.pt -o tmp_out_alignment tmp_out_alignment/junction.bed

#Step 3: output the cleaned alignment file
splam clean -P -o tmp_out_alignment -@ 5
