SAMPLE=SRR1352129_chr9_sub
BAM_FILE=./data/SRR1352129_chr9_sub.bam
MODEL=$3

# Step 1
echo ./Stepls_1_generate_bed.sh $SAMPLE $BAM_FILE
./Step_1_generate_bed.sh $SAMPLE $BAM_FILE

echo python Step_2_get_donors_acceptors_coords.py $SAMPLE
python Step_2_get_donors_acceptors_coords.py $SAMPLE

echo ./Step_3_extract_sequences.sh $SAMPLE
./Step_3_extract_sequences.sh $SAMPLE

echo python Step_4_create_x_y.py $SAMPLE
python Step_4_create_x_y.py $SAMPLE

echo python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/junction.score.bed -m ../src/MODEL/SPLAM_v10/splam_14_scripted.pt
python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/junction.score.bed -m ../src/MODEL/SPLAM_v10/splam_14_scripted.pt

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