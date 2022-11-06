###############################################################################
'''This parser takes as input the .h5 file produced by create_datafile.py and
outputs a .h5 file with datapoints of the form (X, Y), which can be understood
by Keras models.'''
###############################################################################

import h5py
import numpy as np
import sys
import time
from utils_SpliceNN import *
from constants import *

start_time = time.time()

assert sys.argv[1] in ['train', 'test', 'all']

h5f = h5py.File(data_dir + 'datafile'
                + '_' + sys.argv[1]
                + '.h5', 'r')

SEQ = h5f['SEQ'][:]
STRAND = h5f['STRAND'][:]
TX_START = h5f['TX_START'][:]
TX_END = h5f['TX_END'][:]
JN_START = h5f['JN_START'][:]
JN_END = h5f['JN_END'][:]
h5f.close()

SEQ = SEQ.astype('str')
STRAND = STRAND.astype('str')
TX_START = TX_START.astype('str')
TX_END = TX_END.astype('str')
JN_START = JN_START.astype('str')
JN_END = JN_END.astype('str')

h5f2 = h5py.File(data_dir + 'dataset'
                + '_' + sys.argv[1]
                + '.h5', 'w')

CHUNK_SIZE = 100

for i in range(SEQ.shape[0]//CHUNK_SIZE):
    # Each dataset has CHUNK_SIZE genes
    
    if (i+1) == SEQ.shape[0]//CHUNK_SIZE:
        NEW_CHUNK_SIZE = CHUNK_SIZE + SEQ.shape[0]%CHUNK_SIZE
    else:
        NEW_CHUNK_SIZE = CHUNK_SIZE

    X_batch = []
    Y_batch = [[] for t in range(1)]

    for j in range(NEW_CHUNK_SIZE):

        idx = i*CHUNK_SIZE + j

        X, Y = create_datapoints(SEQ[idx], STRAND[idx],
                                 TX_START[idx], TX_END[idx],
                                 JN_START[idx], JN_END[idx])
        print("X.shape: ", X.shape)
        print("Y.shape: ", Y[0].shape)
        X_batch.extend(X)
        print(len(X_batch))
        # print("X_batch.shape: ", np.array(X_batch).shape)
        for t in range(1):
            Y_batch[t].extend(Y[t])

        # print("X_batch.shape: ", X_batch[0].shape)
    X_batch = np.asarray(X_batch).astype('int8')
    for t in range(1):
        Y_batch[t] = np.asarray(Y_batch[t]).astype('int8')

    print("X_batch: ", X_batch.shape)
    print("Y_batch: ", Y_batch[0].shape)

    h5f2.create_dataset('X' + str(i), data=X_batch)
    h5f2.create_dataset('Y' + str(i), data=Y_batch)

h5f2.close()

print(("--- %s seconds ---" % (time.time() - start_time)))

###############################################################################         
