SAMPLE=$1
MODEL=$2

# Step 1
echo Step_1_generate_bed.sh $SAMPLE
./Step_1_generate_bed.sh $SAMPLE

echo python Step_2_get_donors_acceptors_coords.py $SAMPLE ORIGINAL
python Step_2_get_donors_acceptors_coords.py $SAMPLE ORIGINAL

echo ./Step_3_extract_sequences.sh $SAMPLE ORIGINAL
./Step_3_extract_sequences.sh $SAMPLE ORIGINAL

echo python Step_4_create_x_y.py $SAMPLE
python Step_4_create_x_y.py $SAMPLE

echo python script/splam.py -f ../results/800bp/$SAMPLE/INPUTS/input.fa -o ../results/800bp/$SAMPLE/score.bed -m ../src/MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24.pt

python script/splam.py -f ../results/800bp/$SAMPLE/INPUTS/input.fa -o ../results/800bp/$SAMPLE/score.bed -m ../src/MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v32/SpliceNN_24.pt








echo python Step_5_remove_juncs.py $SAMPLE $MODEL
python Step_5_remove_juncs.py $SAMPLE $MODEL

# echo ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE $MODEL
# # ./Step_6_stringtie_assembly.sh $SAMPLE BEFORE $MODEL

# echo ./Step_6_stringtie_assembly.sh $SAMPLE AFTER $MODEL
# # ./Step_6_stringtie_assembly.sh $SAMPLE AFTER $MODEL


# echo ./Step_7_gffcompare.sh $SAMPLE BEFORE  $MODEL ORIGINAL
# # ./Step_7_gffcompare.sh $SAMPLE BEFORE $MODEL ORIGINAL

# echo ./Step_7_gffcompare.sh $SAMPLE AFTER $MODEL ORIGINAL
# # ./Step_7_gffcompare.sh $SAMPLE AFTER $MODEL ORIGINAL