# Step 1: extract introns in the annotation
splam extract lncs_and_exons.gff -o tmp_out_annotation_lncRNA -F feature_transcript.txt

# Step 2: score introns in the annotation
splam score -G ~/Documents/Projects/ref_genome/homo_sapiens/NCBI_Refseq_chr_fixed/GCF_000001405.40_GRCh38.p14_genomic.fna -m ../model/splam_script.pt -o tmp_out_annotation_lncRNA tmp_out_annotation_lncRNA/junction.bed

#Step 3: output statistics of each transcript
splam clean -o tmp_out_annotation_lncRNA -t 0.8
