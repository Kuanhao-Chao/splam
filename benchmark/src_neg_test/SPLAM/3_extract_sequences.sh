SEQ_LEN=800

for target in "GRCm39" "Mmul_10" "NHGRI_mPanTro3" "TAIR10"
#for target in "TAIR10"
do
    echo Generating sequence for $target
    mkdir -p ./3_output/${SEQ_LEN}bp/${target}/

    bedtools getfasta -s -fi ../SPLAM_python/extraction/primates/${target}_genomic.fa -bed ./2_output/${SEQ_LEN}bp/${target}/donor.bed -fo ./3_output/${SEQ_LEN}bp/${target}/donor_seq.fa

    bedtools getfasta -s -fi ../SPLAM_python/extraction/primates/${target}_genomic.fa -bed ./2_output/${SEQ_LEN}bp/${target}/acceptor.bed -fo ./3_output/${SEQ_LEN}bp/${target}/acceptor_seq.fa
done