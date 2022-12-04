if [ $2 = "BEFORE" ]
then
    gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/Dataset/GCF_000001405.40_GRCh38.p14_genomic.primary.gff ../../results/$1/OUTPUT/$3/BEFORE/$1.gtf -o ../../results/$1/OUTPUT/$3/BEFORE/ORIGINAL
fi


if [ $2 = "AFTER" ]
then
    gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SpliceNN/Dataset/GCF_000001405.40_GRCh38.p14_genomic.primary.gff ../../results/$1/OUTPUT/$3/AFTER/$1.gtf -o ../../results/$1/OUTPUT/$3/AFTER/CLEANED
fi
