###############################################################################
'''This code has functions which process the information in the .h5 files
datafile_{}_{}.h5 and convert them into a format usable by Keras.'''
###############################################################################

import numpy as np
import re
from math import ceil
from sklearn.metrics import average_precision_score
from constants import *

assert CL_max % 2 == 0

IN_MAP = np.asarray([[0, 0, 0, 0],
                     [1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
# One-hot encoding of the inputs: 0 is for padding, and 1, 2, 3, 4 correspond
# to A, C, G, T respectively.

OUT_MAP = np.asarray([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1],
                      [0, 0, 0]])
# One-hot encoding of the outputs: 0 is for no splice, 1 is for acceptor,
# 2 is for donor and -1 is for padding.


def ceil_div(x, y):

    return int(ceil(float(x)/y))


def create_datapoints(seq, strand, tx_start, tx_end, jn_start, jn_end):
    # This function first converts the sequence into an integer array, where
    # A, C, G, T, N are mapped to 1, 2, 3, 4, 0 respectively. If the strand is
    # negative, then reverse complementing is done. The splice junctions 
    # are also converted into an array of integers, where 0, 1, 2, -1 
    # correspond to no splicing, acceptor, donor and missing information
    # respectively. It then calls reformat_data and one_hot_encode
    # and returns X, Y which can be used by Keras models.

    main_seq=seq[CL_max//2:-CL_max//2]
    seq = 'N'*(CL_max//2) + main_seq + 'N'*(CL_max//2)
    # Context being provided on the RNA and not the DNA

    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')

    tx_start = int(tx_start)
    tx_end = int(tx_end) 

    # jn_start = [x for x in str(jn_start)]
    # print("jn_start: ", jn_start)
    # print("jn_end  : ", jn_end)
    # print(jn_start[0].split(','))
    jn_start = [x.split(',')[:-1] for x in jn_start]
    jn_end = [x.split(',')[:-1] for x in jn_end]

    jn_start = np.array(jn_start, dtype="int_")
    jn_end = np.array(jn_end, dtype="int_")
    # print("jn_start: ", jn_start)


    if strand == '+':

        X0 = np.asarray(list(map(int, list(seq))))
        Y0 = [-np.ones(tx_end-tx_start+1) for t in range(1)]

        for t in range(1):
            
            if len(jn_start[t]) > 0:
                Y0[t] = np.zeros(tx_end-tx_start+1)
                for c in jn_start[t]:
                    if tx_start <= c <= tx_end:
                        Y0[t][c-tx_start] = 2
                for c in jn_end[t]:
                    if tx_start <= c <= tx_end:
                        Y0[t][c-tx_start] = 1
                    # Ignoring junctions outside annotated tx start/end sites
                     
    elif strand == '-':

        X0 = (5-np.asarray(list(map(int, list(seq[::-1]))))) % 5  # Reverse complement
        Y0 = [-np.ones(tx_end-tx_start+1) for t in range(1)]

        for t in range(1):

            if len(jn_start[t]) > 0:
                Y0[t] = np.zeros(tx_end-tx_start+1)
                for c in jn_end[t]:
                    if tx_start <= c <= tx_end:
                        Y0[t][tx_end-c] = 2
                for c in jn_start[t]:
                    if tx_start <= c <= tx_end:
                        Y0[t][tx_end-c] = 1

    Xd, Yd = reformat_data(X0, Y0)
    X, Y = one_hot_encode(Xd, Yd)

    # X0 = np.array([X0])
    # Y0 = np.array([Y0])
    # X, Y = one_hot_encode(X0, Y0)

    # print(X.shape)
    # print(Y[0].shape)

    return X, Y


def reformat_data(X0, Y0):
    # This function converts X0, Y0 of the create_datapoints function into
    # blocks such that the data is broken down into data points where the
    # input is a sequence of length SL+CL_max corresponding to SL nucleotides
    # of interest and CL_max context nucleotides, the output is a sequence of
    # length SL corresponding to the splicing information of the nucleotides
    # of interest. The CL_max context nucleotides are such that they are
    # CL_max/2 on either side of the SL nucleotides of interest.
    
    num_points = ceil_div(len(Y0[0]), SL)
    # print("num_points: ", num_points)

    Xd = np.zeros((num_points, SL+CL_max))
    Yd = [-np.ones((num_points, SL)) for t in range(1)]

    # print("X0.shape: ", X0.shape)
    # print("Y0.shape: ", Y0[0].shape)

    # print("Xd.shape: ", Xd.shape)
    # print("Yd.shape: ", Yd[0].shape)

    X0 = np.pad(X0, [0, SL], 'constant', constant_values=0)
    Y0 = [np.pad(Y0[t], [0, SL], 'constant', constant_values=-1)
         for t in range(1)]

    # print("X0.shape: ", X0.shape)
    # print("Y0.shape: ", Y0[0].shape)

    for i in range(num_points):
        Xd[i] = X0[SL*i:CL_max+SL*(i+1)]

    for t in range(1):
        for i in range(num_points):
            Yd[t][i] = Y0[t][SL*i:SL*(i+1)]

    return Xd, Yd

def reformat_data_no_chunk(X0, Y0):
    # This function converts X0, Y0 of the create_datapoints function into
    # blocks such that the data is broken down into data points where the
    # input is a sequence of length SL+CL_max corresponding to SL nucleotides
    # of interest and CL_max context nucleotides, the output is a sequence of
    # length SL corresponding to the splicing information of the nucleotides
    # of interest. The CL_max context nucleotides are such that they are
    # CL_max/2 on either side of the SL nucleotides of interest.
    
    num_points = 1
    # # print("num_points: ", num_points)

    Xd = np.zeros((num_points, X0.shape[0]))
    Yd = [-np.ones((num_points, Y0[0].shape[0])) for t in range(1)]

    # print("X0.shape: ", X0.shape)
    # print("Y0.shape: ", Y0[0].shape)

    # # print("Xd.shape: ", Xd.shape)
    # # print("Yd.shape: ", Yd[0].shape)

    # X0 = np.pad(X0, [0, SL], 'constant', constant_values=0)
    # Y0 = [np.pad(Y0[t], [0, SL], 'constant', constant_values=-1)
    #      for t in range(1)]

    # # print("X0.shape: ", X0.shape)
    # # print("Y0.shape: ", Y0[0].shape)

    for i in range(num_points):
        Xd[i] = X0[0:X0.shape[0]]

    for t in range(1):
        for i in range(num_points):
            Yd[t][i] = Y0[t][:Y0[0].shape[0]]

    return Xd, Yd

def clip_datapoints_spliceAI(X, Y, CL, N_GPUS):
    # This function is necessary to make sure of the following:
    # (i) Each time model_m.fit is called, the number of datapoints is a
    # multiple of N_GPUS. Failure to ensure this often results in crashes.
    # (ii) If the required context length is less than CL_max, then
    # appropriate clipping is done below.
    # Additionally, Y is also converted to a list (the .h5 files store
    # them as an array).

    rem = X.shape[0]%N_GPUS
    clip = (CL_max-CL)//2

    if rem != 0 and clip != 0:
        return X[:-rem, clip:-clip], [Y[t][:-rem] for t in range(1)]
    elif rem == 0 and clip != 0:
        return X[:, clip:-clip], [Y[t] for t in range(1)]
    elif rem != 0 and clip == 0:
        return X[:-rem], [Y[t][:-rem] for t in range(1)]
    else:
        return X, [Y[t] for t in range(1)]


def clip_datapoints(X, Y, CL, N_GPUS):
    # This function is necessary to make sure of the following:
    # (i) Each time model_m.fit is called, the number of datapoints is a
    # multiple of N_GPUS. Failure to ensure this often results in crashes.
    # (ii) If the required context length is less than CL_max, then
    # appropriate clipping is done below.
    # Additionally, Y is also converted to a list (the .h5 files store 
    # them as an array).

    rem = X.shape[0]%N_GPUS
    # clip = (CL_max-CL)//2
    # print("rem: ", rem)
    # print("clip: ", CL)

    return X[:, CL:-CL], [Y[t] for t in range(1)]

    # if rem != 0 and clip != 0:
    #     return X[:-rem, clip:-clip], [Y[t][:-rem] for t in range(1)]
    # elif rem == 0 and clip != 0:
    #     return X[:, clip:-clip], [Y[t] for t in range(1)]
    # elif rem != 0 and clip == 0:
    #     return X[:-rem], [Y[t][:-rem] for t in range(1)]
    # else:
    #     return X, [Y[t] for t in range(1)]


def one_hot_encode(Xd, Yd):

    return IN_MAP[Xd.astype('int8')], \
           [OUT_MAP[Yd[t].astype('int8')] for t in range(1)]


def print_topl_statistics(y_true, y_pred):
    # Prints the following information: top-kL statistics for k=0.5,1,2,4,
    # auprc, thresholds for k=0.5,1,2,4, number of true splice sites.
    # print ("y_true: ", y_true.shape)
    # print ("y_pred: ", y_pred.shape)

    idx_true = np.nonzero(y_true == 1)[0]
    # print(("idx_true: ", idx_true))
    argsorted_y_pred = np.argsort(y_pred)
    # print(("argsorted_y_pred: ", argsorted_y_pred))
    # print(("argsorted_y_pred.shape: ", argsorted_y_pred.shape))
    sorted_y_pred = np.sort(y_pred)
    # print(("sorted_y_pred: ", sorted_y_pred))
    # print(("sorted_y_pred.shape: ", sorted_y_pred.shape))

    topkl_accuracy = []
    threshold = []

    for top_length in [0.5, 1, 2, 4, 10]:

        idx_pred = argsorted_y_pred[-int(top_length*len(idx_true)):]
        # print(("idx_pred: ", idx_pred))
        
        # print(("np.size(np.intersect1d(idx_true, idx_pred)): ", np.size(np.intersect1d(idx_true, idx_pred))))
        # print(("float(min(len(idx_pred), len(idx_true))): ", float(min(len(idx_pred), len(idx_true)))))
        topkl_accuracy += [np.size(np.intersect1d(idx_true, idx_pred)) \
                  / float(min(len(idx_pred), len(idx_true)))]
        threshold += [sorted_y_pred[-int(top_length*len(idx_true))]] 

    auprc = average_precision_score(y_true, y_pred)

    # print((("%.4f\t\033[91m%.4f\t\033[0m%.4f\t%.4f\t\033%.4f\t\033[94m%.4f\t\033[0m"
    #       + "\t%d") % (
    #       topkl_accuracy[0], topkl_accuracy[1], topkl_accuracy[2],
    #       topkl_accuracy[3], topkl_accuracy[4], auprc,  len(idx_true))))

    print((("\n%.4f\t\033[91m%.4f\t\033[0m%.4f\t%.4f\t\033%.4f\t\033[94m%.4f\t\033[0m"
          + "%.4f\t%.4f\t%.4f\t%.4f\t%d\n\n") % (
          topkl_accuracy[0], topkl_accuracy[1], topkl_accuracy[2],
          topkl_accuracy[3], topkl_accuracy[4], auprc, threshold[0], threshold[1],
          threshold[2], threshold[3], len(idx_true))))


def print_top_1_statistics(y_true, y_pred):
    
    idx_true = np.nonzero(y_true == 1)[0]
    argsorted_y_pred = np.argsort(y_pred)
    sorted_y_pred = np.sort(y_pred)

    # topkl_accuracy = 

    top_length = 1

    idx_pred = argsorted_y_pred[-int(top_length*len(idx_true)):]
    # print(("idx_true: ", idx_true))
    # print(("idx_pred: ", idx_pred))
    # print(("np.size(np.intersect1d(idx_true, idx_pred)): ", np.size(np.intersect1d(idx_true, idx_pred))))
    # print(("float(min(len(idx_pred), len(idx_true))): ", float(min(len(idx_pred), len(idx_true)))))
    if (len(idx_true) == 0):
        topkl_accuracy = 0
    else:
        topkl_accuracy = np.size(np.intersect1d(idx_true, idx_pred)) \
                    / float(min(len(idx_pred), len(idx_true)))                
    auprc = average_precision_score(y_true, y_pred)
    return topkl_accuracy, auprc

def print_threshold_statistics(y_true, y_pred, threshold, TOTAL_TP, TOTAL_FN, TOTAL_FP):
    idx_true = np.nonzero(y_true == 1)[0]
    idx_pred = np.nonzero(y_pred > threshold)[0]
    # print("idx_true: ", idx_true)
    # print("idx_pred: ", idx_pred)

    LCL_TOTAL_TP = np.size(np.intersect1d(idx_true, idx_pred))
    LCL_TOTAL_FN = len(idx_true) - LCL_TOTAL_TP
    LCL_TOTAL_FP = len(idx_pred) - LCL_TOTAL_TP

    # print("LCL_TOTAL_TP: ", LCL_TOTAL_TP)
    # print("LCL_TOTAL_FN: ", LCL_TOTAL_FN)
    # print("LCL_TOTAL_FP: ", LCL_TOTAL_FP)

    TOTAL_TP += LCL_TOTAL_TP  
    TOTAL_FN += LCL_TOTAL_FN  
    TOTAL_FP += LCL_TOTAL_FP  

    # precision = np.size(np.intersect1d(idx_true, idx_pred)) \
    #             / float(len(idx_pred))        
    # recall = np.size(np.intersect1d(idx_true, idx_pred)) \
    #             / float(len(idx_true)) 

    # print("precision: ", precision)
    # print("recall:    ", recall)
    return TOTAL_TP, TOTAL_FN, TOTAL_FP