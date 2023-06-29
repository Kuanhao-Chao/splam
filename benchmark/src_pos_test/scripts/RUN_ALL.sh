SAMPLE=data
DATABASE=TAIR10

# Step 1
echo "\n\n>> Step 1"
echo Step_1_get_spliceai_donors_acceptors_coords.py $SAMPLE $DATABASE
python Step_1_get_spliceai_donors_acceptors_coords.py $SAMPLE $DATABASE

# Step 2
echo "\n\n>> Step 2"
echo ./Step_3_extract_sequences.sh $SAMPLE $DATABASE
./Step_3_extract_sequences.sh $SAMPLE $DATABASE

# Step 3
echo "\n\n>> Step 3"
echo python Step_4_concatenate_spliceai_input.py $SAMPLE $DATABASE
python Step_4_concatenate_spliceai_input.py $SAMPLE $DATABASE

# Step 4 (technically 6)
# do this manually...