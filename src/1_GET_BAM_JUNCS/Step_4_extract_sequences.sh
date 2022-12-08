THRESHOLD=100
SEQ_LEN=1000

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_junctions/${SEQ_LEN}bp/${THRESHOLD}_juncs/donor.bed -fo ../BAM_junctions/${SEQ_LEN}bp/${THRESHOLD}_juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_junctions/${SEQ_LEN}bp/${THRESHOLD}_juncs/acceptor.bed -fo ../BAM_junctions/${SEQ_LEN}bp/${THRESHOLD}_juncs/acceptor_seq.fa