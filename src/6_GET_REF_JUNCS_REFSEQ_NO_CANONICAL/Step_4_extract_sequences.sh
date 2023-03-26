SEQ_LEN=800
mkdir ./REF_junctions/${SEQ_LEN}bp

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./REF_junctions/${SEQ_LEN}bp/100_juncs/donor.bed -fo ./REF_junctions/${SEQ_LEN}bp/100_juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./REF_junctions/${SEQ_LEN}bp/100_juncs/acceptor.bed -fo ./REF_junctions/${SEQ_LEN}bp/100_juncs/acceptor_seq.fa