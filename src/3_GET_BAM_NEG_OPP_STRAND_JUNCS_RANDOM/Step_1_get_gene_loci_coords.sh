mkdir NEG_rev_junctions

gtf2bed --gtf ../../Dataset/MANE.GRCh38.v1.0.ensembl_genomic.gtf --bed ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.bed

sort -k 1,1 -k2,2n -k3,3n ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.sort.bed

bedtools merge -i ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.sort.bed -s -c 6 -o distinct > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.bed

awk '{if ($4 == "-") {print}}' ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.bed

awk '{if ($4 == "+") {print}}' ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.bed

bedtools intersect -v -a ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.bed -b ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.only.bed

bedtools intersect -v -a ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.bed -b ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.only.bed

cat ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.pos.only.bed ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.neg.only.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.merge.bed

bedtools sort -i ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.merge.bed > ./NEG_rev_junctions/MANE.GRCh38.v1.0.ensembl_genomic.merge.merge.sort.bed