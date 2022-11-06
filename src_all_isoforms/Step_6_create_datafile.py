###############################################################################
'''This parser takes as input the text files canonical_dataset.txt and 
canonical_sequence.txt, and produces a .h5 file datafile_{}_{}.h5,
which will be later processed to create dataset_{}_{}.h5. The file
dataset_{}_{}.h5 will have datapoints of the form (X,Y), and can be
understood by Keras models.'''
###############################################################################

import numpy as np
import re
import sys
import time
import h5py
from constants import *
import random
import math

start_time = time.time()
dt = h5py.string_dtype(encoding='utf-8')

# assert sys.argv[1] in ['train', 'test', 'all']
# assert sys.argv[2] in ['0', '1', 'all']

# if sys.argv[1] == 'train':
#     CHROM_GROUP = ['chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
#                    'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
#                    'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']
# elif sys.argv[1] == 'test':
#     CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9']
# else:
#     CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9',
#                    'chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
#                    'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
#                    'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']

###############################################################################

NAME = []      # Gene symbol
PARALOG = []   # 0 if no paralogs exist, 1 otherwise
CHROM = []     # Chromosome number
STRAND = []    # Strand in which the gene lies (+ or -)
TX_START = []  # Position where transcription starts
TX_END = []    # Position where transcription ends
JN_START = []  # Positions where canonical exons end
JN_END = []    # Positions where canonical exons start
SEQ = []       # Nucleotide sequence

fpr2 = open(sequence, 'r')

with open(splice_table, 'r') as fpr1:
    for line1 in fpr1:

        line2 = fpr2.readline()
        # print("line2: ", line2)

        data1 = re.split('\n|\t', line1)[:-1]
        data2 = re.split('\n|\t|:|-', line2)[:-1]
        # print("data1: ", data1)
        # print("data1[4]: ", data1[4])
        # print("data1[5]: ", data1[5])
        # print("data1[2]: ", data1[2])
        # print("data2[0]: ", data2[0])

        if (data1[2] != data2[0] or int(data1[4]) != int(data2[1])+CL_max//2+1 or int(data1[5]) != int(data2[2])-CL_max//2):
            continue
        assert data1[2] == data2[0]
        assert int(data1[4]) == int(data2[1])+CL_max//2+1
        assert int(data1[5]) == int(data2[2])-CL_max//2

        # if (data1[2] not in CHROM_GROUP):
        #     continue

        # if (sys.argv[2] != data1[1]) and (sys.argv[2] != 'all'):
        #     continue

        NAME.append(data1[0])
        PARALOG.append(int(data1[1]))
        CHROM.append(data1[2])
        STRAND.append(data1[3])
        TX_START.append(data1[4])
        TX_END.append(data1[5])
        JN_START.append(data1[6::2])
        JN_END.append(data1[7::2])
        SEQ.append(data2[3])

fpr1.close()
fpr2.close()

NAME = np.asarray(NAME, dtype=dt)
PARALOG =np.asarray(PARALOG, dtype=dt)
CHROM = np.asarray(CHROM, dtype=dt)
STRAND = np.asarray(STRAND, dtype=dt)
TX_START = np.asarray(TX_START, dtype=dt)
TX_END = np.asarray(TX_END, dtype=dt)
JN_START = np.asarray(JN_START, dtype=dt)
JN_END = np.asarray(JN_END, dtype=dt)
SEQ = np.asarray(SEQ, dtype=dt)

sample_num = len(NAME)


idx_ls = [*range(0,sample_num)]
# print(idx_ls)

random.shuffle(idx_ls)
# print(idx_ls)

TRAIN_RATIO = 0.9
TEST_RATIO = 0.1

train_num = math.ceil(sample_num*0.9)
test_num = sample_num - train_num

idx_train_ls = idx_ls[0:train_num]
idx_test_ls = idx_ls[train_num:]
print("sample_num: ", sample_num)
print(len(idx_train_ls))
print(len(idx_test_ls))

###############################################################################

h5f = h5py.File(data_dir + 'datafile' 
                + '_train'
                + '.h5', 'w')

ds = h5f.create_dataset('NAME', data=np.asarray(NAME[idx_train_ls], dtype=dt) , dtype=dt)

# h5f.create_dataset('PARALOG', data=np.asarray(PARALOG[idx_train_ls]))
h5f.create_dataset('CHROM', data=np.asarray(CHROM[idx_train_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('STRAND', data=np.asarray(STRAND[idx_train_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('TX_START', data=np.asarray(TX_START[idx_train_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('TX_END', data=np.asarray(TX_END[idx_train_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('JN_START', data=np.asarray(JN_START[idx_train_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('JN_END', data=np.asarray(JN_END[idx_train_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('SEQ', data=np.asarray(SEQ[idx_train_ls], dtype=dt) , dtype=dt)

h5f.close()

h5f = h5py.File(data_dir + 'datafile' 
                + '_test'
                + '.h5', 'w')

dt = h5py.string_dtype(encoding='utf-8')
ds = h5f.create_dataset('NAME', data=np.asarray(NAME[idx_test_ls], dtype=dt) , dtype=dt)

# h5f.create_dataset('PARALOG', data=np.asarray(PARALOG[idx_test_ls]))
h5f.create_dataset('CHROM', data=np.asarray(CHROM[idx_test_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('STRAND', data=np.asarray(STRAND[idx_test_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('TX_START', data=np.asarray(TX_START[idx_test_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('TX_END', data=np.asarray(TX_END[idx_test_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('JN_START', data=np.asarray(JN_START[idx_test_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('JN_END', data=np.asarray(JN_END[idx_test_ls], dtype=dt) , dtype=dt)
h5f.create_dataset('SEQ', data=np.asarray(SEQ[idx_test_ls], dtype=dt) , dtype=dt)

h5f.close()



print("--- %s seconds ---" % (time.time() - start_time))

###############################################################################

