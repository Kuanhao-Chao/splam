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
from utils import *
# from utils_SpliceNN import *
# from constants import *

def split_seq_name(seq):
    return seq[1:]

def main(argv):
    start_time = time.time()

    # assert sys.argv[1] in ['train', 'test', 'all']
    # assert sys.argv[2] in ['0', '1', 'all']

    SEQ_LEN = "600"
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

    SEQ = ""
    CHRNAME = ""
    STRAND = ""
    TX_START = ""
    TX_END = ""
    JN_START = ""
    JN_END = ""

    h5f2 = h5py.File("./INPUTS/dataset.h5", 'w')

    pos_faf = "../../src/INPUTS/"+SEQ_LEN+"bp/input_pos/test_pos.shuffle.fa"
    nge_1_faf = "../../src/INPUTS/"+SEQ_LEN+"bp/input_neg_1/test_neg_1.shuffle.fa"
    neg_can_faf = "../../src/INPUTS/"+SEQ_LEN+"bp/input_neg_can/test_neg_can.shuffle.fa"
    neg_noncan_faf = "../../src/INPUTS/"+SEQ_LEN+"bp/input_neg_noncan/test_neg_noncan.shuffle.fa"


    # fa_file = "../../results/spliceAI/"+SEQ_LEN+"bp/"+argv[0]+"/INPUTS/input.fa"

    CONSTANT_SIZE = 23294
    CONSTANT_SIZE_NEG = math.ceil(CONSTANT_SIZE*2/3)


    X_batch = []
    Y_batch = [[] for t in range(1)]
        
    #################################
    ## Processing 'POSITIVE' samples
    #################################
    pidx = 0
    with open(pos_faf, "r") as f:
        print("Processing ", pos_faf)
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
                X_batch.append(X)
                # print(len(X_batch))
                # print("X_batch.shape: ", np.array(X_batch).shape)
                for t in range(1):
                    # Y_batch[t].extend(Y[t])
                    Y_batch[t].append(Y[t])                

            pidx += 1
            if pidx %10000 == 0:
                print("pidx: ", pidx)
            if pidx > CONSTANT_SIZE:
                break

    #################################
    ## Processing 'NEGATIVE_1' samples
    #################################
    n1idx = 0
    with open(nge_1_faf, "r") as f:
        print("Processing ", nge_1_faf)
        lines = f.read().splitlines()
        seq_name = ""
        seq = ""
        for line in lines:
            # print(line)
            if n1idx % 2 == 0:
                seq_name = split_seq_name(line)
            elif n1idx % 2 == 1:
                seq = line
                # print(seq)
                X, Y = create_datapoints(seq, '-')
                X_batch.append(X)
                # print(len(X_batch))
                # print("X_batch.shape: ", np.array(X_batch).shape)
                for t in range(1):
                    # Y_batch[t].extend(Y[t])
                    Y_batch[t].append(Y[t])                

            n1idx += 1
            if n1idx %10000 == 0:
                print("n1idx: ", n1idx)
            if n1idx > CONSTANT_SIZE_NEG:
                break

    #################################
    ## Processing 'NEGATIVE' samples
    #################################
    nidx = 0
    with open(neg_can_faf, "r") as f:
        print("Processing ", neg_can_faf)
        lines = f.read().splitlines()
        seq_name = ""
        seq = ""
        for line in lines:
            # print(line)
            if nidx % 2 == 0:
                seq_name = split_seq_name(line)
            elif nidx % 2 == 1:
                seq = line
                # print(seq)
                X, Y = create_datapoints(seq, '-')
                X_batch.append(X)
                # print(len(X_batch))
                # print("X_batch.shape: ", np.array(X_batch).shape)
                for t in range(1):
                    # Y_batch[t].extend(Y[t])
                    Y_batch[t].append(Y[t])                       
            nidx += 1
            if nidx %10000 == 0:
                print("nidx: ", nidx)
            if nidx > CONSTANT_SIZE_NEG:
                break

    #####################x############
    ## Processing 'Non-canonical NEGATIVE' samples
    #################################
    nnidx = 0
    with open(neg_noncan_faf, "r") as f:
        print("Processing ", neg_noncan_faf)
        lines = f.read().splitlines()
        seq_name = ""
        seq = ""
        for line in lines:
            # print(line)
            if nnidx % 2 == 0:
                seq_name = split_seq_name(line)
            elif nnidx % 2 == 1:
                seq = line
                # print(seq)
                X, Y = create_datapoints(seq, '-')
                X_batch.append(X)
                # print(len(X_batch))
                # print("X_batch.shape: ", np.array(X_batch).shape)
                for t in range(1):
                    # Y_batch[t].extend(Y[t])
                    Y_batch[t].append(Y[t])   

            nnidx += 1
            if nnidx %10000 == 0:
                print("nnidx: ", nnidx)
            if nnidx > CONSTANT_SIZE_NEG:
                break
    # random.shuffle(self.data)

    X_batch = np.asarray(X_batch).astype('int8')
    for t in range(1):
        Y_batch[t] = np.asarray(Y_batch[t]).astype('int8')

    print("\tX_batch: ", X_batch.shape)
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