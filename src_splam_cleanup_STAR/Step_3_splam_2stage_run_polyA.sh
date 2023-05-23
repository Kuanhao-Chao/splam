for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo ../SPLAM_C++/build/splam clean -m ../MODEL/splam_14_script.pt -r ../Dataset/GCF_000001405.40_GRCh38.p14_genomic.fa -o ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/  --paired-removal --2-stage-run ../Dataset/polyA_STAR/$sample/$sample.bamAligned.sortedByCoord.out.bam &
    ../SPLAM_C++/build/splam clean -m ../MODEL/splam_14_script.pt -r ../Dataset/GCF_000001405.40_GRCh38.p14_genomic.fa -o ../results/polyA_STAR/$sample.bamAligned.sortedByCoord.out/  --paired-removal --2-stage-run ../Dataset/polyA_STAR/$sample/$sample.bamAligned.sortedByCoord.out.bam &
done
