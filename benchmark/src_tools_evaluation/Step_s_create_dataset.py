###############################################################################
'''This parser takes as input the .h5 file produced by create_datafile.py and
outputs a .h5 file with datapoints of the form (X, Y), which can be understood
by Keras models.'''
###############################################################################

import h5py
import numpy as np
import sys
import time
import math
import os, sys
from utils import *
# from utils_SpliceNN import *
# from constants import *

def split_seq_name(seq):
    return seq[1:]

def main(argv):
    start_time = time.time()

    # assert sys.argv[1] in ['train', 'test', 'all']
    # assert sys.argv[2] in ['0', '1', 'all']

    type = "test"
    # h5f = h5py.File(data_dir + 'datafile'
    #                 + '_' + sys.argv[1] + '_' + sys.argv[2]
    #                 + '.h5', 'r')

    # SEQ = h5f['SEQ'][:]
    # STRAND = h5f['STRAND'][:]
    # TX_START = h5f['TX_START'][:]
    # TX_END = h5f['TX_END'][:]
    # JN_START = h5f['JN_START'][:]
    # JN_END = h5f['JN_END'][:]
    # h5f.close()
    CONSTANT_SIZE = 10

    SEQ = ""
    CHRNAME = ""
    STRAND = ""
    TX_START = ""
    TX_END = ""
    JN_START = ""
    JN_END = ""

    os.makedirs("./INPUT/", exist_ok=True)
    h5f2 = h5py.File("./INPUT/dataset.h5", 'w')

    da_faf = "./output/spliceai.juncs.seq.fa"

    X_batch = []
    Y_batch = [[] for t in range(1)]
        
    #################################
    ## Processing 'POSITIVE' samples
    #################################
    pidx = 0
    with open(da_faf, "r") as f:
        print("Processing ", da_faf)
        lines = f.read().splitlines()
        seq_name = ""
        seq = ""
        for line in lines:
            # print(line)
            if pidx % 2 == 0:
                seq_name = split_seq_name(line)
            elif pidx % 2 == 1:
                seq = line
                # print(seq)
                X, Y = create_datapoints(seq, '+')
                X_batch.extend(X)
                # print("X_batch: ", X_batch)
                print(len(X_batch))
                print("X_batch.shape: ", np.array(X_batch).shape)
                for t in range(1):
                    # Y_batch[t].extend(Y[t])
                    Y_batch[t].extend(Y[t])
                    # print("Y_batch[t].shape: ", np.array(Y_batch[t]).shape)

            pidx += 1
            if pidx %10000 == 0:
                print("pidx: ", pidx)
            # if pidx > 100:
            #     break
    # print(X_batch)
    # print("np.asarray(X_batch): ", np.array(X_batch).shape)
    X_batch = np.asarray(X_batch).astype('int8')
    print("X_batch: ", X_batch.shape)
    for t in range(1):
        Y_batch[t] = np.asarray(Y_batch[t]).astype('int8')
        print("Y_batch[t]: ", Y_batch[t].shape)

    # print("\tX_batch: ", X_batch.shape)
    # print("\tY_batch: ", Y_batch[0].shape)

    h5f2.create_dataset('X', data=X_batch)
    h5f2.create_dataset('Y', data=Y_batch)
    # CHUNK_IDX += 1
    # X_batch = []
    # Y_batch = [[] for t in range(1)]
    h5f2.close()

    print(("--- %s seconds ---" % (time.time() - start_time)))


if __name__ == "__main__":
    main(sys.argv[1:])