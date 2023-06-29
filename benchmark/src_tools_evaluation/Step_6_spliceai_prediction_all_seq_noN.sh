TARGET=$1
SPLICEAI_VERSION=$2
for (( c=1000; c<=22000; c+=1000 ))
do
    echo ">> python Step_6_spliceai_prediction_all_seq.py $c noN $TARGET $SPLICEAI_VERSION"
    python Step_6_spliceai_prediction_all_seq.py $c noN $TARGET $SPLICEAI_VERSION
    echo "\n\n"
done
