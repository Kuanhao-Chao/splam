
# bedtools getfasta -s -fi ../Dataset/GRCh38_p14_genomic_ucsc.fa -bed $1/juncs/donor.bed -fo $1/juncs/donor_seq.fa

# bedtools getfasta -s -fi ../Dataset/GRCh38_p14_genomic_ucsc.fa -bed $1/juncs/acceptor.bed -fo $1/juncs/acceptor_seq.fa

bedtools getfasta -s -fi ../Dataset/${2}_genomic.fa -bed $1/juncs/donor.bed -fo $1/juncs/donor_seq.fa

bedtools getfasta -s -fi ../Dataset/${2}_genomic.fa -bed $1/juncs/acceptor.bed -fo $1/juncs/acceptor_seq.fa