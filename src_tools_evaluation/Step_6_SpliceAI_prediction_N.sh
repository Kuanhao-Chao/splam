TARGET=$1
for (( c=200; c<=20400; c+=200 ))
# for (( c=20; c<=100; c+=20 ))
do
    echo python Step_6_SpliceAI_prediction.py $c N $TARGET
    python Step_6_SpliceAI_prediction.py $c N $TARGET
done