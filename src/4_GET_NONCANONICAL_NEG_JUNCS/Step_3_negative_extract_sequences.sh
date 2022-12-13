SEQ_LENGTH=600

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../NEG_noncan_junctions/${SEQ_LENGTH}bp/donor.bed -fo ../NEG_noncan_junctions/${SEQ_LENGTH}bp/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../NEG_noncan_junctions/${SEQ_LENGTH}bp/acceptor.bed -fo ../acceptor_seq.fa
NEG_noncan_junctions/${SEQ_LENGTH}bp/acceptor_seq.fa