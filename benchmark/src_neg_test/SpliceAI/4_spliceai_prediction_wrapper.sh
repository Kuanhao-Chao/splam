SPLICEAI_VERSION=$1

for target in "GRCm39" "Mmul_10" "NHGRI_mPanTro3" "TAIR10"
do 
    echo ./4_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    ./4_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    echo \n\n

    echo ./4_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    ./4_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    echo \n\n
done