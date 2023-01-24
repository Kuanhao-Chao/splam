SAMPLE=$1
# TARGET="NEG_noncan_junctions"

for OUTPUTFILE in "./OUTPUT/pos/" "./OUTPUT/neg_can/" "./OUTPUT/neg_noncan/" "./OUTPUT/neg_1/"
do  
    echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/spliceai.juncs.bed -fo $OUTPUTFILE/spliceai.juncs.seq.fa"

    bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/spliceai.juncs.bed -fo $OUTPUTFILE/spliceai.juncs.seq.fa

    echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/splam.juncs.donor.bed -fo $OUTPUTFILE/splam.juncs.donor.seq.fa"

    bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/splam.juncs.donor.bed -fo $OUTPUTFILE/splam.juncs.donor.seq.fa

    echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/splam.juncs.acceptor.bed -fo $OUTPUTFILE/splam.juncs.acceptor.seq.fa"

    bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed $OUTPUTFILE/splam.juncs.acceptor.bed -fo $OUTPUTFILE/splam.juncs.acceptor.seq.fa
done