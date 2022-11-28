# TARGET="NEG_junctions"
TARGET="NEG_noncan_junctions"

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../$TARGET/donor.bed -fo ../$TARGET/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../$TARGET/acceptor.bed -fo ../$TARGET/acceptor_seq.fa