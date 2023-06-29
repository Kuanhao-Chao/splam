DATABASE=TAIR10

# Step 1
echo "\n\n>> Step 1"
echo python 1_get_spliceai_donors_acceptors_coords.py $DATABASE
python 1_get_spliceai_donors_acceptors_coords.py $DATABASE

# Step 2
echo "\n\n>> Step 2"
echo ./2_extract_sequences.sh $DATABASE
./2_extract_sequences.sh $DATABASE

# Step 3
echo "\n\n>> Step 3"
echo python 3_concatenate_spliceai_input.py $SAMPLE $DATABASE
python 3_concatenate_spliceai_input.py $SAMPLE $DATABASE

# Step 4 (technically 6)
# do this manually...