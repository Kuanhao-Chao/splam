gtf2bed --gtf ../../Dataset/gencode.v43.chr_patch_hapl_scaff.annotation.gtf --bed ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.bed

sort -k 1,1 -k2,2n -k3,3n ../../Dataset/gencode.v43.chr_patch_hapl_scaff.annotation.bed > ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.sort.bed

bedtools merge -i ../../Dataset/gencode.v43.chr_patch_hapl_scaff.annotation.sort.bed -s -c 6 -o distinct > ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.bed

awk '{if ($4 == "-") {print}}' ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.bed > ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.neg.bed

awk '{if ($4 == "+") {print}}' ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.bed > ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.pos.bed

bedtools intersect -v -a ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.pos.bed -b ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.neg.bed > ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.pos.only.bed

bedtools intersect -v -a ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.neg.bed -b ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.pos.bed > ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.neg.only.bed

bedtools intersect -wa -wb -a ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.pos.only.bed -b BAM_junctions/junctions_1.bed -sorted -filenames > ./BAM_junctions/neg_hits.bed

# bedtools intersect -wa -wb -a ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.pos.only.bed -b BAM_junctions/junctions_1.bed -sorted -filenames > ./BAM_junctions/neg_hits.bed

bedtools intersect -wa -wb -a ./BAM_junctions/gencode.v43.chr_patch_hapl_scaff.annotation.merge.neg.only.bed -b BAM_junctions/junctions_1.bed -sorted -filenames > ./BAM_junctions/pos_hits.bed