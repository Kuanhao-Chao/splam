SAMPLE=tests
DATABASE=refseq

# Step 1
#echo ">> Step 1"
#echo "../SPLAM_C++/build/splam j-extract -o ./$SAMPLE $BAM_FILE"
#../SPLAM_C++/build/splam j-extract -o ./$SAMPLE $BAM_FILE

# Step 2
echo "\n\n>> Step 2"
echo "python Step_2_get_donors_acceptors_coords.py $SAMPLE ${DATABASE}_parsed.bed"
python Step_2_get_donors_acceptors_coords.py $SAMPLE ${DATABASE}_parsed.bed

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
echo python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/${DATABASE}.score.bed -m ../MODEL/splam_14_script.pt
python splam_predict.py -f $SAMPLE/INPUTS/input.fa -o $SAMPLE/${DATABASE}.score.bed -m ../MODEL/splam_14_script.pt