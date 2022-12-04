bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../../results/$1/juncs/donor.bed -fo ../../results/$1/juncs/donor_seq.fa

bedtools getfasta -s -fi ../../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ../../results/$1/juncs/acceptor.bed -fo ../../results/$1/juncs/acceptor_seq.fa