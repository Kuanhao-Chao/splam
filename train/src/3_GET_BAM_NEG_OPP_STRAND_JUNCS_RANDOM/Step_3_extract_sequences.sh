threshold=5
SEQ_LEN=800
bedtools getfasta -s -fi ../../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./NEG_rev_junctions/${SEQ_LEN}bp/donor/donor.bed -fo ./NEG_rev_junctions/${SEQ_LEN}bp/donor/donor_seq.fa

bedtools getfasta -s -fi ../../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./NEG_rev_junctions/${SEQ_LEN}bp/acceptor/acceptor.bed -fo ./NEG_rev_junctions/${SEQ_LEN}bp/acceptor/acceptor_seq.fa