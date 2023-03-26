gtf2bed --gtf ../../Dataset/MANE.GRCh38.v1.0.ensembl_genomic.gtf --bed ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.bed

sort -k 1,1 -k2,2n -k3,3n ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.bed > ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.sort.bed

bedtools merge -i ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.sort.bed -s -c 6 -o distinct > ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.bed

awk '{if ($4 == "-") {print}}' ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.bed > ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.bed

awk '{if ($4 == "+") {print}}' ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.bed > ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.bed

bedtools intersect -v -a ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.bed -b ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.bed > ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.only.bed

bedtools intersect -v -a ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.bed -b ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.bed > ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.only.bed

bedtools intersect -wa -wb -a ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.only.bed -b BAM_junctions/junctions_1.bed -sorted -filenames > ./BAM_junctions/neg_hits.bed

bedtools intersect -wa -wb -a ./BAM_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.only.bed -b BAM_junctions/junctions_1.bed -sorted -filenames > ./BAM_junctions/pos_hits.bed