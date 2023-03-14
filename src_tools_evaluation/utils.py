import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from math import ceil
import numpy as np
from sklearn.metrics import average_precision_score

JUNC_START = 150
JUNC_END = 450
CL_MAX = 10000

def ceil_div(x, y):

    return int(ceil(float(x)/y))

#######################################
# This is for Conformer model 
#######################################
def create_datapoints(seq, strand):
    # print("seq length: ", len(seq))
    # seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0').replace('K', '0').replace('R', '0')
    jn_start = 0
    jn_end = len(seq)-10000-1
    # print("Donor: ", seq[CL_MAX//2+jn_start: CL_MAX//2+jn_start+2])
    # print("Donor: ", seq[CL_MAX//2+jn_end-2: CL_MAX//2+jn_end])

    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = [np.zeros(len(seq)-CL_MAX) for t in range(1)]

    if strand == '+':
        for t in range(1):        
            Y0[t][jn_start] = 2
            Y0[t][jn_end] = 1
            # print(" Y0[t][jn_end] : ",  Y0[t].shape)

    # X0, Y0 = reformat_data(X0, Y0)
    X, Y = one_hot_encode(X0, Y0)

    # print("X.shape: ", X.shape)
    # print("Y.shape: ", Y[0].shape)

    return X, Y

def create_datapoints_spliceai_N(seq, strand):
    # print("seq length: ", len(seq))
    # seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0').replace('K', '0').replace('R', '0')
    jn_start = 200
    jn_end = len(seq)-10000-1-200
    # print("Donor: ", seq[CL_MAX//2+jn_start: CL_MAX//2+jn_start+2])
    # print("Donor: ", seq[CL_MAX//2+jn_end-2: CL_MAX//2+jn_end])

    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = [np.zeros(len(seq)-CL_MAX) for t in range(1)]

    if strand == '+':
        for t in range(1):        
            Y0[t][jn_start] = 2
            Y0[t][jn_end] = 1
            # print(" Y0[t][jn_end] : ",  Y0[t].shape)

    # X0, Y0 = reformat_data(X0, Y0)
    X, Y = one_hot_encode(X0, Y0)

    # print("X.shape: ", X.shape)
    # print("Y.shape: ", Y[0].shape)

    return X, Y



def reformat_data(X0, Y0):
    # This function converts X0, Y0 of the create_datapoints function into
    # blocks such that the data is broken down into data points where the
    # input is a sequence of length SL+CL_max corresponding to SL nucleotides
    # of interest and CL_max context nucleotides, the output is a sequence of
    # length SL corresponding to the splicing information of the nucleotides
    # of interest. The CL_max context nucleotides are such that they are
    # CL_max/2 on either side of the SL nucleotides of interest.

    SL = 5000
    CL_max = 10000
    num_points = ceil_div(len(Y0[0]), SL)
    # print("num_points: ", num_points)

    Xd = np.zeros((num_points, SL+CL_max))
    Yd = [-np.ones((num_points, SL)) for t in range(1)]

    print("X0.shape: ", X0.shape)
    print("Y0.shape: ", Y0[0].shape)

    # print("Xd.shape: ", Xd.shape)
    # print("Yd.shape: ", Yd[0].shape)

    X0 = np.pad(X0, [0, SL], 'constant', constant_values=0)
    Y0 = [np.pad(Y0[t], [0, SL], 'constant', constant_values=-1)
         for t in range(1)]

    # print("X0.shape: ", X0.shape)
    # print("Y0.shape: ", Y0[0].shape)

    for i in range(num_points):
        # print("SL*i:CL_max+SL*(i+1) = ", SL*i, " - ", CL_max+SL*(i+1))
        Xd[i] = X0[SL*i:CL_max+SL*(i+1)]

    for t in range(1):
        for i in range(num_points):
            Yd[t][i] = Y0[t][SL*i:SL*(i+1)]

    print("final Xd.shape: ", Xd.shape)
    print("final Yd.shape: ", Yd[0].shape)

    return Xd, Yd




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

