SPLICEAI_VERSION=$1

for target in "GRCm39" "Mmul_10" "NHGRI_mPanTro3" "TAIR10"
do 
    echo ./Step_6_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    ./Step_6_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    echo \n\n

    echo ./Step_6_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    ./Step_6_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    echo \n\n
done