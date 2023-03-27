SAMPLE=$1
MODEL=$2

# Step 1
# echo Step_1_generate_bed.sh $SAMPLE
# ./Step_1_generate_bed.sh $SAMPLE

echo python Step_2_get_donors_acceptors_coords.py $SAMPLE ORIGINAL
python Step_2_get_donors_acceptors_coords.py $SAMPLE ORIGINAL

echo ./Step_3_extract_sequences.sh $SAMPLE ORIGINAL
./Step_3_extract_sequences.sh $SAMPLE ORIGINAL

echo python Step_4_create_x_y.py $SAMPLE
python Step_4_create_x_y.py $SAMPLE

echo python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/score.bed -m ../src/MODEL/SPLAM_v10/splam_14_scripted.pt
python splam_predict.py  -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/score.bed -m ../src/MODEL/SPLAM_v10/splam_14_scripted.pt

# echo python Step_5_remove_juncs.py $SAMPLE $MODEL
# python Step_5_remove_juncs.py $SAMPLE $MODEL

# echo ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE $MODEL
# # ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE $MODEL

# echo ./Step_6_stringtie_assembly.sh $SAMPLE AFTER $MODEL
# # ./Step_6_stringtie_assembly.sh $SAMPLE AFTER $MODEL


# echo ./Step_7_gffcompare.sh $SAMPLE BEFORE  $MODEL ORIGINAL
# # ./Step_7_gffcompare.sh $SAMPLE BEFORE $MODEL ORIGINAL

# echo ./Step_7_gffcompare.sh $SAMPLE AFTER $MODEL ORIGINAL
# # ./Step_7_gffcompare.sh $SAMPLE AFTER $MODEL ORIGINAL