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

start_time = time.time()

assert sys.argv[1] in ['train', 'test', 'all']
assert sys.argv[2] in ['0', '1', 'all']

if sys.argv[1] == 'train':
    CHROM_GROUP = ['chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
                   'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
                   'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']
elif sys.argv[1] == 'test':
    CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9']
else:
    CHROM_GROUP = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9',
                   'chr11', 'chr13', 'chr15', 'chr17', 'chr19', 'chr21',
                   'chr2', 'chr4', 'chr6', 'chr8', 'chr10', 'chr12',
                   'chr14', 'chr16', 'chr18', 'chr20', 'chr22', 'chrX', 'chrY']

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

print(fpr2)
with open(splice_table, 'r') as fpr1:
    for line1 in fpr1:

        line2 = fpr2.readline()

        data1 = re.split('\n|\t', line1)[:-1]
        data2 = re.split('\n|\t|:|-', line2)[:-1]

        assert data1[2] == data2[0]
        assert int(data1[4]) == int(data2[1])+CL_max//2+1
        assert int(data1[5]) == int(data2[2])-CL_max//2

        if (data1[2] not in CHROM_GROUP):
            continue

        if (sys.argv[2] != data1[1]) and (sys.argv[2] != 'all'):
            continue

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

###############################################################################

h5f = h5py.File(data_dir + 'datafile' 
                + '_' + sys.argv[1] + '_' + sys.argv[2]
                + '.h5', 'w')

dt = h5py.string_dtype(encoding='utf-8')
ds = h5f.create_dataset('NAME', data=np.asarray(NAME, dtype=dt) , dtype=dt)

h5f.create_dataset('PARALOG', data=np.asarray(PARALOG))
h5f.create_dataset('CHROM', data=np.asarray(CHROM, dtype=dt) , dtype=dt)
h5f.create_dataset('STRAND', data=np.asarray(STRAND, dtype=dt) , dtype=dt)
h5f.create_dataset('TX_START', data=np.asarray(TX_START, dtype=dt) , dtype=dt)
h5f.create_dataset('TX_END', data=np.asarray(TX_END, dtype=dt) , dtype=dt)
h5f.create_dataset('JN_START', data=np.asarray(JN_START, dtype=dt) , dtype=dt)
h5f.create_dataset('JN_END', data=np.asarray(JN_END, dtype=dt) , dtype=dt)
h5f.create_dataset('SEQ', data=np.asarray(SEQ, dtype=dt) , dtype=dt)

h5f.close()

print("--- %s seconds ---" % (time.time() - start_time))

###############################################################################

