TARGET=$1
for (( c=1000; c<=10000; c+=1000 ))
do
    echo ">> python Step_6_spliceai_prediction_all_seq.py $c N $TARGET"
    python Step_6_spliceai_prediction_all_seq.py $c N $TARGET
    echo "\n\n"
done