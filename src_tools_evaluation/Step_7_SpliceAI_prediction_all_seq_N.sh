TARGET=$1
for (( c=200; c<=20200; c+=200 ))
# for (( c=20; c<=100; c+=20 ))
do
    echo python Step_7_SpliceAI_prediction_all_seq.py $c N $TARGET
    python Step_7_SpliceAI_prediction_all_seq.py $c N $TARGET
done