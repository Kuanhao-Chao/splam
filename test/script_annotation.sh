# Step 1: extract introns in the annotation
#splam extract refseq_110_GRCh38_chr_fixed.gff -o tmp_out_annotation -F feature_gene.txt
splam extract MANE.GRCh38.v1.1.subset.gff -o tmp_out_annotation_MANE -F feature_gene.txt

# Step 2: score introns in the annotation
splam score -G chr9_subset.fa -m ../model/splam_static.pt -o tmp_out_annotation_MANE tmp_out_annotation_MANE/junction.bed

#Step 3: output statistics of each transcript
splam clean -o tmp_out_annotation_MANE -t 0.8
