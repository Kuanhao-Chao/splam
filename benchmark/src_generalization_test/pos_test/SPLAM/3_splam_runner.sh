for DATABASE in "GRCm39" "Mmul_10" "NHGRI_mPanTro3" "TAIR10"
do
    echo Generating scores for $DATABASE
    mkdir -p ./3_output/${DATABASE}/

    echo python 3_splam_predict.py -f ./2_output/${DATABASE}/input_neg_random.fa -o ./3_output/${DATABASE}/${DATABASE}.score.bed -m ../../../model/splam_script.pt
    python 3_splam_predict.py -f ./2_output/${DATABASE}/input_neg_random.fa -o ./3_output/${DATABASE}/${DATABASE}.score.bed -m ../../../model/splam_script.pt
done
