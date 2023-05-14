First, update the constants.py file:
- ref_genome: path of the genome.fa file (hg19/GRCh37)

Then, use the following commands:

./Step_1_grab_sequence.sh

python Step_2_create_datafile.py train all
python Step_2_create_datafile.py test 0

python Step_3_create_dataset.py train all
python Step_3_create_dataset.py test 0

python Step_4_verify_h5_file.py
python Step_5_test_model.py 10000 test_1



qsub script_train.sh 80 1
qsub script_train.sh 80 2
qsub script_train.sh 80 3
qsub script_train.sh 80 4
qsub script_train.sh 80 5

qsub script_test.sh 10000

# The code was tested using keras==2.0.5 and tensorflow==1.4.1
