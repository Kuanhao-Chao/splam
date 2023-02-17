SAMPLE=$1
MODEL=$2

mkdir ../results/800bp/$SAMPLE/assembly
touch ../results/800bp/$MODEL.txt

echo samtools sort ../results/800bp/$1/cleaned.bam -o ../results/800bp/$1/cleaned.sort.bam
samtools sort ../results/800bp/$1/cleaned.bam -o ../results/800bp/$1/cleaned.sort.bam

echo ./Step_1_stringtie_assembly.sh $SAMPLE BEFORE
./Step_1_stringtie_assembly.sh $SAMPLE BEFORE

echo ./Step_1_stringtie_assembly.sh $SAMPLE AFTER
./Step_1_stringtie_assembly.sh $SAMPLE AFTER


echo ./Step_2_gffcompare.sh $SAMPLE BEFORE HISAT2
./Step_2_gffcompare.sh $SAMPLE BEFORE HISAT2

echo ./Step_2_gffcompare.sh $SAMPLE AFTER HISAT2
./Step_2_gffcompare.sh $SAMPLE AFTER HISAT2