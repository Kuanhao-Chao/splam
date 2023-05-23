# polyA
for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do

    mkdir ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/
    mkdir ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/
    mkdir ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/BEFORE
    mkdir ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/AFTER

    gffcompare -r ../Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff  ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/assembly/BEFORE/$sample.gtf -o ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/BEFORE/res &
    
    gffcompare -r ../Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff  ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/$sample.gtf -o ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/AFTER/res &
done

# ribozero
for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287; do

    mkdir ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/
    mkdir ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/
    mkdir ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/BEFORE
    mkdir ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/AFTER

    gffcompare -r ../Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff  ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/assembly/BEFORE/$sample.gtf -o ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/BEFORE/res &

    gffcompare -r ../Dataset/refseq_GCF_000001405.40_GRCh38.p14_genomic.gff  ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/assembly/AFTER/$sample.gtf -o ../results/ribozero_STAR/$sample.bamAligned.sortedByCoord.out/gffcompare/refseq_ucsc/AFTER/res &
done