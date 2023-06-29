TARGET=$1
SPLICEAI_VERSION=$2
for (( c=1000; c<=50000; c+=1000 ))
do
    echo ">> python 4_spliceai_prediction_all_seq.py $c noN $TARGET $SPLICEAI_VERSION"
    python 4_spliceai_prediction_all_seq.py $c noN $TARGET $SPLICEAI_VERSION
    echo \n\n
done
