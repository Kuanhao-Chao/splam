# polyA
for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do
    echo mkdir ../results/polyA/$sample/intron_matcher
    mkdir ../results/polyA/$sample/intron_matcher

    echo mkdir ../results/polyA/$sample/intron_matcher/BEFORE/
    mkdir ../results/polyA/$sample/intron_matcher/BEFORE/

    echo ../../Intron-Matcher/build/intron_matcher -G ../Dataset/hg38.ncbiRefSeq.gtf -o ../results/polyA/$sample/intron_matcher/BEFORE/res.txt ../Dataset/polyA/$sample/$sample.bam &
    ../../Intron-Matcher/build/intron_matcher -G ../Dataset/hg38.ncbiRefSeq.gtf -o ../results/polyA/$sample/intron_matcher/BEFORE/res.txt ../Dataset/polyA/$sample/$sample.bam &
done

# ribozero
for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287; do
    echo mkdir ../results/ribozero/$sample/intron_matcher
    mkdir ../results/ribozero/$sample/intron_matcher

    echo mkdir ../results/ribozero/$sample/intron_matcher/BEFORE/
    mkdir ../results/ribozero/$sample/intron_matcher/BEFORE/

    echo ../../Intron-Matcher/build/intron_matcher -G ../Dataset/hg38.ncbiRefSeq.gtf -o ../results/ribozero/$sample/intron_matcher/BEFORE/res.txt ../Dataset/ribozero/$sample/$sample.bam &
    ../../Intron-Matcher/build/intron_matcher -G ../Dataset/hg38.ncbiRefSeq.gtf -o ../results/ribozero/$sample/intron_matcher/BEFORE/res.txt ../Dataset/ribozero/$sample/$sample.bam &
done
