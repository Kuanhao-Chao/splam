TARGET=$1
SPLICEAI_VERSION=$2
for (( c=500; c<=25000; c+=500 ))
do
    echo python 5_spliceai_prediction_all_seq.py $c N $TARGET $SPLICEAI_VERSION
    python 5_spliceai_prediction_all_seq.py $c N $TARGET $SPLICEAI_VERSION
    echo -e "\n"
done
 