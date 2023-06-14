SPLICEAI_VERSION=$1
for target in "pos_MANE" "pos_ALTS" "neg_1" "neg_random"
#for target in "pos_MANE" "pos_ALTS"
do
    echo ./Step_6_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    ./Step_6_spliceai_prediction_all_seq_noN.sh $target $SPLICEAI_VERSION &
    echo "\n\n"

    echo ./Step_6_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    ./Step_6_spliceai_prediction_all_seq_N.sh $target $SPLICEAI_VERSION &
    echo "\n\n"
done
