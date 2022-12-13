<<<<<<< HEAD
SEQ_LEN=600
mkdir ../BAM_REF_Intersection/${SEQ_LEN}bp
=======
SEQ_LEN=1000
>>>>>>> 1189cf671af213485edd35714556970d3b41c338

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_REF_Intersection/${SEQ_LEN}bp/100_juncs/donor.bed -fo ../BAM_REF_Intersection/${SEQ_LEN}bp/100_juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_REF_Intersection/${SEQ_LEN}bp/100_juncs/acceptor.bed -fo ../BAM_REF_Intersection/${SEQ_LEN}bp/100_juncs/acceptor_seq.fa