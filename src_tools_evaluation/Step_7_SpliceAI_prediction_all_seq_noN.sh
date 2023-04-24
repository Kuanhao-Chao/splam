TARGET=$1
for (( c=1000; c<=25000; c+=1000 ))
do
    echo ">> python Step_7_SpliceAI_prediction_all_seq.py $c noN $TARGET"
    # python Step_7_SpliceAI_prediction_all_seq.py $c noN $TARGET
    echo "\n\n"
done