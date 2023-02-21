SAMPLE=$1
TARGET=$2

SEQ_LEN=800
if [ $TARGET = "BEFORE" ]
then
    mkdir ./results/$SAMPLE/assembly
    mkdir ./results/$SAMPLE/assembly/BEFORE

    echo "samtools index ../Dataset/$SAMPLE/$SAMPLE.bam"
    samtools index ../Dataset/$SAMPLE/$SAMPLE.bam
    
    echo "stringtie -o ./results/$SAMPLE/assembly/BEFORE/$SAMPLE.gtf ../Dataset/$SAMPLE/$SAMPLE.bam"
    stringtie -o ./results/$SAMPLE/assembly/BEFORE/$SAMPLE.gtf ../Dataset/$SAMPLE/$SAMPLE.bam
fi

if [ $TARGET = "AFTER" ]
then
    mkdir ./results/$SAMPLE/assembly
    mkdir ./results/$SAMPLE/assembly/AFTER

    echo "samtools index /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/$SAMPLE/cleaned.fix.sort.bam"
    samtools index /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/$SAMPLE/cleaned.fix.sort.bam
    
    echo "stringtie -o ./results/$SAMPLE/assembly/AFTER/$SAMPLE.gtf /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/$SAMPLE/cleaned.fix.sort.bam"
    stringtie -o ./results/$SAMPLE/assembly/AFTER/$SAMPLE.gtf /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/$SAMPLE/cleaned.fix.sort.bam
fi