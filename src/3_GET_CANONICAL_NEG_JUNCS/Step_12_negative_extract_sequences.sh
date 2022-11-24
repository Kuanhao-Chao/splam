bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./NEG_junctions/donor.bed -fo ./NEG_junctions/donor_seq.fa

bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./NEG_junctions/acceptor.bed -fo ./NEG_junctions/acceptor_seq.fa