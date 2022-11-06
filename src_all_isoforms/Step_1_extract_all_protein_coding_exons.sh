awk '{if($3=="exon"){print}}' ../Dataset/GCF_000001405.40_GRCh38.p14_genomic.gtf | grep "transcript_biotype\ \"mRNA\"" > ../Dataset/all_protein_coding_exons.gtf
