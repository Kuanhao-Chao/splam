bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_junctions/100_juncs/donor.bed -fo ../BAM_junctions/100_juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_junctions/100_juncs/acceptor.bed -fo ../BAM_junctions/100_juncs/acceptor_seq.fa