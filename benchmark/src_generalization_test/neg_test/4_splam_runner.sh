for DATABASE in "GRCm39" "Mmul_10" "NHGRI_mPanTro3" "TAIR10"
do
    echo Generating scores for $DATABASE
    mkdir -p ./4_output/${DATABASE}/

    echo python 4_splam_predict.py -f ./2_output/${DATABASE}/input_neg_random.fa -o ./4_output/${DATABASE}/score.bed -m ../../../model/splam_script.pt
    python 4_splam_predict.py -f ./2_output/${DATABASE}/input_neg_random.fa -o ./4_output/${DATABASE}/score.bed -m ../../../model/splam_script.pt
done