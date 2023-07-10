DATABASE=$1

echo Generating sequence for $DATABASE
mkdir -p ./2_output/

echo bedtools getfasta -s -fi ../../SPLAM_python/extraction/primates/${DATABASE}_genomic.fa -bed ./1_output/${DATABASE}_spliceai.bed -fo ./2_output/${DATABASE}_spliceai_seq.fa
bedtools getfasta -s -fi ../../SPLAM_python/extraction/primates/${DATABASE}_genomic.fa -bed ./1_output/${DATABASE}_spliceai.bed -fo ./2_output/${DATABASE}_spliceai_seq.fa