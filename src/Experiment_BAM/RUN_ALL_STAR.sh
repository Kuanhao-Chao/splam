# SAMPLE=$1

# Step 1
echo python Step_2_get_donors_acceptors_coords_STAR.py $SAMPLE
python Step_2_get_donors_acceptors_coords_STAR.py $SAMPLE

echo ./Step_3_extract_sequences_STAR.sh $SAMPLE
./Step_3_extract_sequences_STAR.sh $SAMPLE

echo python Step_4_create_x_y.py $SAMPLE
python Step_4_create_x_y.py $SAMPLE

# echo python SpliceNN_prediction.py $SAMPLE
# python SpliceNN_prediction.py $SAMPLE

# echo python Step_5_remove_juncs.py $SAMPLE
# python Step_5_remove_juncs.py $SAMPLE

# echo ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE
# ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE

# echo ./Step_6_stringtie_assembly.sh $SAMPLE AFTER
# ./Step_6_stringtie_assembly.sh $SAMPLE AFTER



# echo ./Step_7_gffcompare_STAR.sh $SAMPLE BEFORE
# ./Step_7_gffcompare_STAR.sh $SAMPLE BEFORE

# echo ./Step_7_gffcompare.sh $SAMPLE AFTER
# ./Step_7_gffcompare_STAR.sh $SAMPLE AFTER