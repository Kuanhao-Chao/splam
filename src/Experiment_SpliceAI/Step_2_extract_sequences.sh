SEQ_LEN=600

if [ $2 = "STAR" ]
then
bedtools getfasta -s -fi ../../Dataset/GRCh38_latest_genomic.fna -bed ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/donor.bed -fo ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/GRCh38_latest_genomic.fna -bed ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/acceptor.bed -fo ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/acceptor_seq.fa

else
bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/donor.bed -fo ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/acceptor.bed -fo ../../results/spliceAI/${SEQ_LEN}bp/$1/juncs/acceptor_seq.fa

fi


