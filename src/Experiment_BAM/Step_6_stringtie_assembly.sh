mkdir ../../results/$1/OUTPUT/$3/BEFORE
mkdir ../../results/$1/OUTPUT/$3/AFTER

if [ $2 = "BEFORE" ]
then
    stringtie -o ../../results/$1/OUTPUT/$3/BEFORE/$1.gtf ../../Dataset/$1/$1.bam
fi

if [ $2 = "AFTER" ]
then
    samtools index ../../results/$1/OUTPUT/$3/BAM/$1.cleaned.bam
    stringtie -o ../../results/$1/OUTPUT/$3/AFTER/$1.gtf ../../results/$1/OUTPUT/$3/BAM/$1.cleaned.bam
fi