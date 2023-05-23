SAMPLE=$1
MODEL=$2

mkdir ./results
mkdir ./results/$SAMPLE/
mkdir ./results/$SAMPLE/assembly
touch ./results/$MODEL.txt

# echo samtools sort /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/$SAMPLE/cleaned.bam -o ../results/800bp/$SAMPLE/cleaned.sort.bam
# samtools sort /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_SPLAM/build/$SAMPLE/cleaned.bam -o ../results/800bp/$SAMPLE/cleaned.sort.bam



#echo ./Step_1_stringtie_assembly.sh $SAMPLE BEFORE
#./Step_1_stringtie_assembly.sh $SAMPLE BEFORE

echo ./Step_1_stringtie_assembly.sh $SAMPLE AFTER
./Step_1_stringtie_assembly.sh $SAMPLE AFTER


#echo ./Step_2_gffcompare.sh $SAMPLE BEFORE HISAT2
#./Step_2_gffcompare.sh $SAMPLE BEFORE HISAT2

echo "./Step_2_gffcompare.sh $SAMPLE AFTER HISAT2"
./Step_2_gffcompare.sh $SAMPLE AFTER HISAT2
