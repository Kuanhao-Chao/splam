# polyA
for sample in R2826 R2835 R2839 R2845 R2855 R2857 R2869 R2874 R2894 R2895 ; do

    mkdir ../results/polyA/$sample/gffcompare/
    mkdir ../results/polyA/$sample/gffcompare/chess/
    mkdir ../results/polyA/$sample/gffcompare/chess/BEFORE
    mkdir ../results/polyA/$sample/gffcompare/chess/AFTER

    gffcompare -r ../Dataset/chess3.0.1.gtf  ../results/polyA/$sample/assembly/BEFORE/$sample.gtf -o ../results/polyA/$sample/gffcompare/chess/BEFORE/res &
    
    gffcompare -r ../Dataset/chess3.0.1.gtf  ../results/polyA/$sample/assembly/AFTER/$sample.gtf -o ../results/polyA/$sample/gffcompare/chess/AFTER/res &
done

# ribozero
for sample in R12258  R12260  R12263  R12265  R12266  R12277  R12278  R12280  R12285  R12287; do

    mkdir ../results/ribozero/$sample/gffcompare/
    mkdir ../results/ribozero/$sample/gffcompare/chess/
    mkdir ../results/ribozero/$sample/gffcompare/chess/BEFORE
    mkdir ../results/ribozero/$sample/gffcompare/chess/AFTER

    gffcompare -r ../Dataset/chess3.0.1.gtf  ../results/ribozero/$sample/assembly/BEFORE/$sample.gtf -o ../results/ribozero/$sample/gffcompare/chess/BEFORE/res &

    gffcompare -r ../Dataset/chess3.0.1.gtf  ../results/ribozero/$sample/assembly/AFTER/$sample.gtf -o ../results/ribozero/$sample/gffcompare/chess/AFTER/res &
done