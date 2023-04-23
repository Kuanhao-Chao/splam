for (( c=1000; c<=13000; c+=1000 ))
# for (( c=20; c<=100; c+=20 ))
do
    echo python Step_7_SpliceAI_prediction_all_seq.py $c N
    python Step_7_SpliceAI_prediction_all_seq.py $c N
    echo "\n\n"
done