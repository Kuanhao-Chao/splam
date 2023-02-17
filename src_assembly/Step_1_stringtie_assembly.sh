mkdir ../results/800bp/$1/assembly/BEFORE
mkdir ../results/800bp/$1/assembly/AFTER

SEQ_LEN=800
if [ $2 = "BEFORE" ]
then
    samtools index ../Dataset/$1/$1.bam
    stringtie -o ../results/800bp/$1/assembly/BEFORE/$1.gtf ../Dataset/$1/$1.bam
fi

if [ $2 = "AFTER" ]
then
    samtools index ../results/800bp/$1/cleaned.sort.bam
    stringtie -o ../results/800bp/$1/assembly/AFTER/$1.gtf ../results/800bp/$1/cleaned.sort.bam
fi