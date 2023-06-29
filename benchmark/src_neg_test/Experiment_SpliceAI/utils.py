import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from sklearn.metrics import average_precision_score

SEQ_LEN = "600"
JUNC_START = 150
JUNC_END = 450
CL_MAX = 10000












# def create_datapoints(seq, strand, tx_start, tx_end, jn_start, jn_end):
#     # This function first converts the sequence into an integer array, where
#     # A, C, G, T, N are mapped to 1, 2, 3, 4, 0 respectively. If the strand is
#     # negative, then reverse complementing is done. The splice junctions 
#     # are also converted into an array of integers, where 0, 1, 2, -1 
#     # correspond to no splicing, acceptor, donor and missing information
#     # respectively. It then calls reformat_data and one_hot_encode
#     # and returns X, Y which can be used by Keras models.

#     main_seq=seq[CL_max//2:-CL_max//2]
#     seq = 'N'*(CL_max//2) + main_seq + 'N'*(CL_max//2)
#     # Context being provided on the RNA and not the DNA

#     seq = seq.upper().replace('A', '1').replace('C', '2')
#     seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')

#     tx_start = int(tx_start)
#     tx_end = int(tx_end) 

#     # jn_start = [x for x in str(jn_start)]
#     # print("jn_start: ", jn_start)
#     # print("jn_end  : ", jn_end)
#     # print(jn_start[0].split(','))
#     jn_start = [x.split(',')[:-1] for x in jn_start]
#     jn_end = [x.split(',')[:-1] for x in jn_end]

#     jn_start = np.array(jn_start, dtype="int_")
#     jn_end = np.array(jn_end, dtype="int_")
#     # print("jn_start: ", jn_start)


#     if strand == '+':

#         X0 = np.asarray(list(map(int, list(seq))))
#         Y0 = [-np.ones(tx_end-tx_start+1) for t in range(1)]

#         for t in range(1):
            
#             if len(jn_start[t]) > 0:
#                 Y0[t] = np.zeros(tx_end-tx_start+1)
#                 for c in jn_start[t]:
#                     if tx_start <= c <= tx_end:
#                         Y0[t][c-tx_start] = 2
#                 for c in jn_end[t]:
#                     if tx_start <= c <= tx_end:
#                         Y0[t][c-tx_start] = 1
#                     # Ignoring junctions outside annotated tx start/end sites
                     
#     elif strand == '-':

#         X0 = (5-np.asarray(list(map(int, list(seq[::-1]))))) % 5  # Reverse complement
#         Y0 = [-np.ones(tx_end-tx_start+1) for t in range(1)]

#         for t in range(1):

#             if len(jn_start[t]) > 0:
#                 Y0[t] = np.zeros(tx_end-tx_start+1)
#                 for c in jn_end[t]:
#                     if tx_start <= c <= tx_end:
#                         Y0[t][tx_end-c] = 2
#                 for c in jn_start[t]:
#                     if tx_start <= c <= tx_end:
#                         Y0[t][tx_end-c] = 1

#     Xd, Yd = reformat_data_no_chunk(X0, Y0)
#     X, Y = one_hot_encode(Xd, Yd)

#     # X0 = np.array([X0])
#     # Y0 = np.array([Y0])
#     # X, Y = one_hot_encode(X0, Y0)

#     # print(X.shape)
#     # print(Y[0].shape)

#     return X, Y




#######################################
# This is for Conformer model 
#######################################
def create_datapoints(seq, strand):
    # print("seq length: ", len(seq))
    seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0').replace('K', '0').replace('R', '0')
    jn_start = JUNC_START
    jn_end = JUNC_END
    # print("Donor: ", seq[CL_MAX//2+jn_start: CL_MAX//2+jn_start+2])
    # print("Donor: ", seq[CL_MAX//2+jn_end-2: CL_MAX//2+jn_end])

    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = [np.zeros(int(SEQ_LEN)) for t in range(1)]

    if strand == '+':
        for t in range(1):        
            Y0[t][jn_start] = 2
            Y0[t][jn_end] = 1
            # print(" Y0[t][jn_end] : ",  Y0[t].shape)
    X, Y = one_hot_encode(X0, Y0)
    return X, Y


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

def one_hot_encode(Xd, Yd):
    return IN_MAP[Xd.astype('int8')], \
           [OUT_MAP[Yd[t].astype('int8')] for t in range(1)]


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

def print_threshold_statistics(y_true, y_pred, threshold, TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN):
    idx_true = np.nonzero(y_true == 1)[0]
    idx_pred = np.nonzero(y_pred > threshold)[0]

    LCL_TOTAL_TP = np.size(np.intersect1d(idx_true, idx_pred))
    LCL_TOTAL_FN = len(idx_true) - LCL_TOTAL_TP
    LCL_TOTAL_FP = len(idx_pred) - LCL_TOTAL_TP
    LCL_TOTAL_TN = len(y_true) - LCL_TOTAL_TP - LCL_TOTAL_FN - LCL_TOTAL_FP

    # print("LCL_TOTAL_TP: ", LCL_TOTAL_TP)
    # print("LCL_TOTAL_FN: ", LCL_TOTAL_FN)
    # print("LCL_TOTAL_FP: ", LCL_TOTAL_FP)

    TOTAL_TP += LCL_TOTAL_TP
    TOTAL_FN += LCL_TOTAL_FN
    TOTAL_FP += LCL_TOTAL_FP
    TOTAL_TN += LCL_TOTAL_TN

    # precision = np.size(np.intersect1d(idx_true, idx_pred)) \
    #             / float(len(idx_pred))        
    # recall = np.size(np.intersect1d(idx_true, idx_pred)) \
    #             / float(len(idx_true)) 

    # print("precision: ", precision)
    # print("recall:    ", recall)
    return TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN, LCL_TOTAL_TP, LCL_TOTAL_FN, LCL_TOTAL_FP, LCL_TOTAL_TN

def print_junc_statistics(D_YL, A_YL, D_YP, A_YP, threshold, TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN):

    label_junc_idx = (D_YL[:, 150]==1) & (A_YL[:, 450]==1)
    label_nonjunc_idx = (D_YL[:, 150]==0) & (A_YL[:, 450]==0)
    predict_junc_idx = (D_YP[:, 150]>=threshold) & (A_YP[:, 450]>=threshold)
    predict_nonjunc_idx = (D_YP[:, 150]<threshold) | (A_YP[:, 450]<threshold)

    idx_true = np.nonzero(label_junc_idx == True)[0]
    idx_pred = np.nonzero(predict_junc_idx == True)[0]

    # print("idx_true: ", idx_true)
    # print("idx_pred: ", idx_pred)

    LCL_TOTAL_TP = np.size(np.intersect1d(idx_true, idx_pred))
    LCL_TOTAL_FN = len(idx_true) - LCL_TOTAL_TP
    LCL_TOTAL_FP = len(idx_pred) - LCL_TOTAL_TP

    # LCL_TOTAL_TN = np.size(np.intersect1d(label_nonjunc_idx, predict_nonjunc_idx))
    LCL_TOTAL_TN = len(D_YL) - LCL_TOTAL_TP - LCL_TOTAL_FN - LCL_TOTAL_FP

    # print("LCL_TOTAL_TP: ", LCL_TOTAL_TP)
    # print("LCL_TOTAL_FN: ", LCL_TOTAL_FN)
    # print("LCL_TOTAL_FP: ", LCL_TOTAL_FP)
    # print("LCL_TOTAL_TN: ", LCL_TOTAL_TN)

    TOTAL_TP += LCL_TOTAL_TP
    TOTAL_FN += LCL_TOTAL_FN
    TOTAL_FP += LCL_TOTAL_FP
    TOTAL_TN += LCL_TOTAL_TN

    # precision = np.size(np.intersect1d(idx_true, idx_pred)) \
    #             / float(len(idx_pred))        
    # recall = np.size(np.intersect1d(idx_true, idx_pred)) \
    #             / float(len(idx_true)) 

    # print("precision: ", precision)
    # print("recall:    ", recall)
    return TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN, LCL_TOTAL_TP, LCL_TOTAL_FN, LCL_TOTAL_FP, LCL_TOTAL_TN