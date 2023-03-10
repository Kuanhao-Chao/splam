# TARGET="NEG_junctions"
TARGET="NEG_noncan_junctions"

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../TEST/$TARGET/donor.bed -fo ../TEST/$TARGET/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../TEST/$TARGET/acceptor.bed -fo ../TEST/$TARGET/acceptor_seq.fa