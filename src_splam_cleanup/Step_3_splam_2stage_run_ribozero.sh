for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287 ; do
    echo ../SPLAM_C++/build/SPLAM clean -m ../MODEL/splam_14_script.pt -r ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -o ../results/ribozero/$sample/  --paired-removal --2-stage-run ../Dataset/ribozero/$sample/$sample.bam &
    ../SPLAM_C++/build/SPLAM clean -m ../MODEL/splam_14_script.pt -r ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -o ../results/ribozero/$sample/  --paired-removal --2-stage-run ../Dataset/ribozero/$sample/$sample.bam &
done
