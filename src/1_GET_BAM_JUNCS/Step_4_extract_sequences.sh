bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_junctions/donor.bed -fo ../BAM_junctions/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../BAM_junctions/acceptor.bed -fo ../BAM_junctions/acceptor_seq.fa