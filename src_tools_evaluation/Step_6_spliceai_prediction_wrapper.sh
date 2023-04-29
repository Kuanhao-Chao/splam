# for target in "pos_MANE" "pos_ALTS" "neg_1" "neg_random"
for target in "pos_MANE" "pos_ALTS"
do
    echo ./Step_6_spliceai_prediction_all_seq_noN.sh $target &
    ./Step_6_spliceai_prediction_all_seq_noN.sh $target &
    echo "\n\n"

    echo ./Step_6_spliceai_prediction_all_seq_N.sh $target &
    ./Step_6_spliceai_prediction_all_seq_N.sh $target &
    echo "\n\n"
done