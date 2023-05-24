SAMPLE=SRR1352129_chr9_sub
BAM_FILE=./data/SRR1352129_chr9_sub.bam
MODEL=$3

# Step 1
#echo ">> Step 1"
#echo "../SPLAM_C++/build/splam j-extract -o ./$SAMPLE $BAM_FILE"
#../SPLAM_C++/build/splam j-extract -o ./$SAMPLE $BAM_FILE

# Step 2
echo "\n\n>> Step 2"
echo "python Step_2_get_donors_acceptors_coords.py $SAMPLE"
python Step_2_get_donors_acceptors_coords.py $SAMPLE

# Step 3
echo "\n\n>> Step 3"
echo ./Step_3_extract_sequences.sh $SAMPLE
./Step_3_extract_sequences.sh $SAMPLE

# Step 4
echo "\n\n>> Step 4"
echo python Step_4_create_x_y.py $SAMPLE
python Step_4_create_x_y.py $SAMPLE

# Step 5
echo "\n\n>> Step 5"
echo python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/junction.score.bed -m ../MODEL/splam_14_script.pt
python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/junction.score.bed -m ../MODEL/splam_14_script.pt

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
