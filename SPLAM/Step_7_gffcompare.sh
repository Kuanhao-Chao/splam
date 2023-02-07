SEQ_LEN=800

if [ $4 = "STAR" ]
then
    if [ $2 = "BEFORE" ]
    then
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/BEFORE/$1.gtf -o ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/BEFORE/ORIGINAL
    fi


    if [ $2 = "AFTER" ]
    then
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/AFTER/$1.gtf -o ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/AFTER/CLEANED
    fi

else

    if [ $2 = "BEFORE" ]
    then
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/BEFORE/$1.gtf -o ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/BEFORE/ORIGINAL
    fi


    if [ $2 = "AFTER" ]
    then
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/AFTER/$1.gtf -o ../../results/${SEQ_LEN}bp/$1/OUTPUT/$3/AFTER/CLEANED
    fi
fi

