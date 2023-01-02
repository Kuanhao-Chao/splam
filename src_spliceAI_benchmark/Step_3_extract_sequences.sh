SAMPLE=$1
# TARGET="NEG_noncan_junctions"

echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./output/spliceai.juncs.bed -fo ./output/spliceai.juncs.seq.fa"

bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./output/spliceai.juncs.bed -fo ./output/spliceai.juncs.seq.fa

echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./output/splam.juncs.donor.bed -fo ./output/splam.juncs.donor.seq.fa"

bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./output/splam.juncs.donor.bed -fo ./output/splam.juncs.donor.seq.fa

echo "bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./output/splam.juncs.acceptor.bed -fo ./output/splam.juncs.acceptor.seq.fa"

bedtools getfasta -s -fi ../Dataset/hg38_p12_ucsc.no_alts.no_fixs.fa -bed ./output/splam.juncs.acceptor.bed -fo ./output/splam.juncs.acceptor.seq.fa