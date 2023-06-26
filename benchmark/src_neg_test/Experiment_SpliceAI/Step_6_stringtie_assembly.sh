mkdir ../../results/spliceAI/$1/OUTPUT/$3/BEFORE
mkdir ../../results/spliceAI/$1/OUTPUT/$3/AFTER

SEQ_LEN=600
if [ $2 = "BEFORE" ]
then
    stringtie -o ../../results/spliceAI/${SEQ_LEN}bp/$1/OUTPUT/$3/BEFORE/$1.gtf ../../Dataset/$1/$1.bam
fi

if [ $2 = "AFTER" ]
then
    samtools index ../../results/spliceAI/${SEQ_LEN}bp/$1/OUTPUT/$3/BAM/$1.cleaned.bam
    stringtie -o ../../results/spliceAI/${SEQ_LEN}bp/$1/OUTPUT/$3/AFTER/$1.gtf ../../results/spliceAI/${SEQ_LEN}bp/$1/OUTPUT/$3/BAM/$1.cleaned.bam
fi