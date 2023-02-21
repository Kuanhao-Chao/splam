SEQ_LEN=800

if [ $3 = "STAR" ]
then
    if [ $2 = "BEFORE" ]
    then
        echo "gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/BEFORE/$1.gtf -o ./results/$1/assembly/BEFORE/BEFORE"
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/BEFORE/$1.gtf -o ./results/$1/assembly/BEFORE/BEFORE
    fi


    if [ $2 = "AFTER" ]
    then
        echo "gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/AFTER/$1.gtf -o ./results/$1/assembly/AFTER/AFTER"
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/AFTER/$1.gtf -o ./results/$1/assembly/AFTER/AFTER
    fi

else

    if [ $2 = "BEFORE" ]
    then

        echo "gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/BEFORE/$1.gtf -o ./results/$1/assembly/BEFORE/BEFORE"
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/BEFORE/$1.gtf -o ./results/$1/assembly/BEFORE/BEFORE
    fi


    if [ $2 = "AFTER" ]
    then
        echo "gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/AFTER/$1.gtf -o ./results/$1/assembly/AFTER/AFTER"
        gffcompare -r /Users/chaokuan-hao/Documents/Projects/PR_SPLAM/Dataset/GRCh38_latest_genomic_chr_replaced.gff  ./results/$1/assembly/AFTER/$1.gtf -o ./results/$1/assembly/AFTER/AFTER
    fi
fi

