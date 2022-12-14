TARGET="NEG_junctions"
SEQ_LENGTH=800
# TARGET="NEG_noncan_junctions"

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../$TARGET/${SEQ_LENGTH}bp/donor.bed -fo ../$TARGET/${SEQ_LENGTH}bp/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../$TARGET/${SEQ_LENGTH}bp/acceptor.bed -fo ../$TARGET/${SEQ_LENGTH}bp/acceptor_seq.fa