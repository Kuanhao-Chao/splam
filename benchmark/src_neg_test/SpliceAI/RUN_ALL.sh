DATABASE=Mmul_10
BASE_LEN=800

# Step 1
echo "\n\n>> Step 1"
echo python 1_get_spliceai_donors_acceptors_coords.py $DATABASE $BASE_LEN
python 1_get_spliceai_donors_acceptors_coords.py $DATABASE $BASE_LEN

# Step 2
echo "\n\n>> Step 2"
echo ./2_extract_sequences.sh $DATABASE
./2_extract_sequences.sh $DATABASE

# Step 3
echo "\n\n>> Step 3"
echo python 3_concatenate_spliceai_input.py $DATABASE
python 3_concatenate_spliceai_input.py $DATABASE