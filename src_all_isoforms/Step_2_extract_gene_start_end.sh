awk '{if($3=="gene"){print}}' ../Dataset/GCF_000001405.40_GRCh38.p14_genomic.gtf | grep "gene_biotype\ \"protein_coding\"" > ../Dataset/all_protein_coding_genes.gtf
