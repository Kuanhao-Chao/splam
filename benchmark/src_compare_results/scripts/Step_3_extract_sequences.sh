SAMPLE=$1
DATABASE=$2

echo bedtools getfasta -s -fi ../Dataset/${DATABASE}_genomic.fa -bed ../output/${SAMPLE}/${DATABASE}_spliceai.bed -fo ../output/${SAMPLE}/${DATABASE}_spliceai_seq.fa
bedtools getfasta -s -fi ../Dataset/${DATABASE}_genomic.fa -bed ../output/${SAMPLE}/${DATABASE}_spliceai.bed -fo ../output/${SAMPLE}/${DATABASE}_spliceai_seq.fa
