SAMPLE=$1
MODEL=$2

# # Step 1
# echo python Step_2_get_donors_acceptors_coords.py $SAMPLE STAR
# python Step_2_get_donors_acceptors_coords.py $SAMPLE STAR

# echo ./Step_3_extract_sequences.sh $SAMPLE STAR
# ./Step_3_extract_sequences.sh $SAMPLE STAR

# echo python Step_4_create_x_y.py $SAMPLE
# python Step_4_create_x_y.py $SAMPLE

# echo python SpliceNN_prediction.py $SAMPLE $MODEL
# python SpliceNN_prediction.py $SAMPLE $MODEL

# echo python Step_5_remove_juncs.py $SAMPLE $MODEL
# python Step_5_remove_juncs.py $SAMPLE $MODEL

# echo ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE $MODEL
# ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE $MODEL

# echo ./Step_6_stringtie_assembly.sh $SAMPLE AFTER $MODEL
# ./Step_6_stringtie_assembly.sh $SAMPLE AFTER $MODEL


echo ./Step_7_gffcompare.sh $SAMPLE BEFORE  $MODEL STAR
./Step_7_gffcompare.sh $SAMPLE BEFORE $MODEL STAR

echo ./Step_7_gffcompare.sh $SAMPLE AFTER $MODEL STAR
./Step_7_gffcompare.sh $SAMPLE AFTER $MODEL STAR