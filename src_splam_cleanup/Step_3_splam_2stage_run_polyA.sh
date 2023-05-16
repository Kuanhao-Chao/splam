for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo ../SPLAM_C++/build/SPLAM clean -m ../MODEL/splam_14_script.pt -r ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -o ../results/polyA/$sample/  --paired-removal --2-stage-run ../Dataset/polyA/$sample/$sample.bam &
    ../SPLAM_C++/build/SPLAM clean -m ../MODEL/splam_14_script.pt -r ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -o ../results/polyA/$sample/  --paired-removal --2-stage-run ../Dataset/polyA/$sample/$sample.bam &
done
