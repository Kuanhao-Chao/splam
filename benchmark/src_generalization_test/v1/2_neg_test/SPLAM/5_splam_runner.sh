DATABASE=NHGRI_mPanTro3
SEQ_LEN=800

echo Generating sequence for $DATABASE
mkdir -p ./5_output/${SEQ_LEN}bp/${DATABASE}/

echo python 5_splam_predict.py -f ./4_output/${SEQ_LEN}bp/${DATABASE}/input_neg_random.fa -o ./5_output/${SEQ_LEN}bp/${DATABASE}/${DATABASE}.score.bed -m ../../model/splam_script.pt
python 5_splam_predict.py -f ./4_output/${SEQ_LEN}bp/${DATABASE}/input_neg_random.fa -o ./5_output/${SEQ_LEN}bp/${DATABASE}/${DATABASE}.score.bed -m ../../model/splam_script.pt