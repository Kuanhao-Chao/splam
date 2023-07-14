target=$1

# for target in "GRCm39" "Mmul_10" "NHGRI_mPanTro3" "TAIR10"
# do 
for SPLICEAI_VERSION in 1 2 3 4 5
do
    echo ./5_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    ./5_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    echo -e "\n"

    echo ./5_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    ./5_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    echo -e "\n"
done
# done