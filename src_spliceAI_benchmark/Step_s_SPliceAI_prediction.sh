for (( c=400; c<=12000; c+=400 ))
# for (( c=20; c<=100; c+=20 ))
do
    echo python Step_6_N_SpliceAI_prediction.py $c
    python Step_6_N_SpliceAI_prediction.py $c
done