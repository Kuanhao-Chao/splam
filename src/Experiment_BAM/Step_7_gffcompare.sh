if [ $2 = "BEFORE" ]
then
    gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/hg38c_protein_and_lncRNA.gtf  ../../results/$1/OUTPUT/$3/BEFORE/$1.gtf -o ../../results/$1/OUTPUT/$3/BEFORE/ORIGINAL
fi


if [ $2 = "AFTER" ]
then
    gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_MultiStringTie/hg38c_protein_and_lncRNA.gtf  ../../results/$1/OUTPUT/$3/AFTER/$1.gtf -o ../../results/$1/OUTPUT/$3/AFTER/CLEANED
fi
