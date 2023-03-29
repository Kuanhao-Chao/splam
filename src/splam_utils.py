###############################################################################
'''This code has functions which process the information in the .h5 files
datafile_{}_{}.h5 and convert them into a format usable by Keras.'''
###############################################################################

import torch
from torch.optim.lr_scheduler import LambdaLR
from torch.optim import Optimizer, AdamW
from torch.nn import CrossEntropyLoss, BCELoss, BatchNorm1d, ModuleList
# from torch.nn import Module, BatchNorm1d, LazyBatchNorm1d, ReLU, LeakyReLU, Conv1d, LazyConv1d, ModuleList, Softmax, Sigmoid, Flatten, Dropout2d, Linear

import numpy as np
import re
import math
from math import ceil
from sklearn.metrics import average_precision_score
from splam_constant import *
from SPLAM import *

SEQ_LEN = 800
# fix random seed
def same_seeds(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

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
           
def one_hot_encode_classifier(Xd):
    return IN_MAP[Xd.astype('int8')]
           
#######################################
# This is for Conformer model 
#######################################
def create_datapoints(seq, strand):
    # seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0').replace('K', '0').replace('R', '0')
    jn_start = JUNC_START
    jn_end = JUNC_END
    # print("Donor: ", seq[CL_MAX//2+jn_start: CL_MAX//2+jn_start+2])
    # print("Donor: ", seq[CL_MAX//2+jn_end-2: CL_MAX//2+jn_end])

    #######################################
    # predicting pb for every bp
    #######################################
    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = [np.zeros(SEQ_LEN) for t in range(1)]
    if strand == '+':
        for t in range(1):        
            Y0[t][jn_start] = 2
            Y0[t][jn_end] = 1
    X, Y = one_hot_encode(X0, Y0)
    return X, Y

    # #######################################
    # # predicting splice / non-splice
    # #######################################
    # X0 = np.asarray(list(map(int, list(seq))))
    # Y0 = 0
    # if strand == '+':
    #     Y0 = 1
    # X = one_hot_encode_classifier(X0)
    # return X, Y0

def get_cosine_schedule_with_warmup(
      optimizer: Optimizer,
      num_warmup_steps: int,
      num_training_steps: int,
      num_cycles: float = 0.5,
      last_epoch: int = -1,
    ):
    """
    Create a schedule with a learning rate that decreases following the values of the cosine function between the
    initial lr set in the optimizer to 0, after a warmup period during which it increases linearly between 0 and the
    initial lr set in the optimizer.

    Args:
    optimizer (:class:`~torch.optim.Optimizer`):
      The optimizer for which to schedule the learning rate.
    num_warmup_steps (:obj:`int`):
      The number of steps for the warmup phase.
    num_training_steps (:obj:`int`):
      The total number of training steps.
    num_cycles (:obj:`float`, `optional`, defaults to 0.5):
      The number of waves in the cosine schedule (the defaults is to just decrease from the max value to 0
      following a half-cosine).
    last_epoch (:obj:`int`, `optional`, defaults to -1):
      The index of the last epoch when resuming training.

    Return:
    :obj:`torch.optim.lr_scheduler.LambdaLR` with the appropriate schedule.
    """

    def lr_lambda(current_step):
        # Warmup
        if current_step < num_warmup_steps:
            return float(current_step) / float(max(1, num_warmup_steps))
        # decadence
        progress = float(current_step - num_warmup_steps) / float(
          max(1, num_training_steps - num_warmup_steps)
        )
        return max(
          0.0, 0.5 * (1.0 + math.cos(math.pi * float(num_cycles) * 2.0 * progress))
        )

    return LambdaLR(optimizer, lr_lambda, last_epoch)

# def model_fn(DNAs, labels, model, criterion):
#     """Forward a batch through the model."""
#     outs = model(DNAs)
#     # print("outs: ", outs.size())
#     # print("labels: ", labels.size())
#     # print("DNAs: ", DNAs.size())

#     # labels = labels.sum(axis=1)
#     print("outs: ", outs)
#     print("labels: ", labels)
#     # print("labels.unsqueeze(1): ", labels.unsqueeze(1))

#     loss = criterion(outs, labels)
#     # print("loss: ", loss)
#     acc = binary_acc(outs, labels)
#     # print("acc: ", acc)

#     return loss, acc, outs

def binary_acc(y_pred, y_test):
    y_pred_tag = torch.round(torch.sigmoid(y_pred))

    # print("y_pred_tag: ", y_pred_tag)

    correct_results_sum = (y_pred_tag == y_test).sum().float()
    acc = correct_results_sum/y_test.shape[0]
    acc = torch.round(acc * 100)
    
    return acc


def get_accuracy(y_prob, y_true):
    assert y_true.ndim == 1 and y_true.size() == y_prob.size()
    y_prob = y_prob > 0.5
    return (y_true == y_prob).sum().item() / y_true.size(0)

# def model_fn_eval(DNAs, labels, model, criterion):
#     """Forward a batch through the model."""
#     outs = model(DNAs)
#     outs = torch.flatten(outs)
#     # print("outs: ", outs.size())
#     # print("outs: ", outs)
#     # print("labels: ", labels.size())
#     # print("DNAs: ", DNAs.size())

#     # labels = labels.sum(axis=1)
#     # print("labels: ", labels.size())

#     loss, accuracy = categorical_crossentropy_2d(labels, outs, criterion)
#     # print("loss    : ", loss)
#     # print("accuracy: ", accuracy)
#     return loss, accuracy, outs

def model_fn(DNAs, labels, model, criterion):
    """Forward a batch through the model."""
    outs = model(DNAs)
    # outs = torch.flatten(outs)
    # print("outs: ", outs.size())
    # print("outs: ", outs)
    # print("labels: ", labels.size())
    # print("DNAs: ", DNAs.size())

    # labels = labels.sum(axis=1)
    # print("labels: ", labels.size())

    loss = categorical_crossentropy_2d(labels, outs, criterion)
    # print("loss    : ", loss)
    # print("accuracy: ", accuracy)
    return loss, outs


def weighted_binary_cross_entropy(output, target, weights=None):
        
    if weights is not None:
        assert len(weights) == 2
        loss = weights[1] * (target * torch.log(output+1e-10)) + \
               weights[0] * ((1 - target) * torch.log(1 - output+1e-10))
    else:
        loss = target * torch.log(output+1e-10) + (1 - target) * torch.log(1 - output+1e-10)

    return torch.neg(torch.mean(loss))

def categorical_crossentropy_2d(y_true, y_pred, criterion):
    # print("y_true: ", y_true)
    # print("y_pred: ", y_pred)
    
    ########################################
    # splice / nonsplice classifier
    ########################################
    # # weights = torch.FloatTensor([1.0, 4.0]) 
    # weights = torch.FloatTensor([1.0, 20.0]) 
    # return weighted_binary_cross_entropy(y_pred, y_true, weights), get_accuracy(y_pred, y_true)
    
    ########################################
    # splice / nonsplice classifier
    ########################################
    SEQ_WEIGHT = 4
    # IMBALANCE_WEIGHT = 1
    # WEIGHT = SEQ_WEIGHT * IMBALANCE_WEIGHT

    # prod = output[:,0]*target
    # return -prod[prod<0].sum()
    # print("y_true: ", y_true)
    # print("y_pred: ", y_pred)
    # print("Loss: ", - torch.mean(y_true[:, 0, :]*torch.log(y_pred[:, 0, :]+1e-10)
    #                     + y_true[:, 1, :]*torch.log(y_pred[:, 1, :]+1e-10)
    #                     + y_true[:, 2, :]*torch.log(y_pred[:, 2, :]+1e-10)))
    # print("y_true[:, 0, :]: ", y_true[:, 0, :])
    # print("y_pred[:, 0, :]: ", y_pred[:, 0, :])
    # print("y_true[:, 1, :]: ", y_true[:, 1, :])
    # print("y_pred[:, 1, :]: ", y_pred[:, 1, :])
    # print("y_true[:, 2, :]: ", y_true[:, 2, :])
    # print("y_pred[:, 2, :]: ", y_pred[:, 2, :])

    # # This is focal loss
    # gamma = 2
    # return - torch.mean(y_true[:, 0, :] * torch.mul( torch.pow( torch.sub(1, y_pred[:, 0, :]), gamma ), torch.log(y_pred[:, 0, :]+1e-10) )
    #                     + SEQ_WEIGHT * y_true[:, 1, :] * torch.mul( torch.pow( torch.sub(1, y_pred[:, 1, :]), gamma ), torch.log(y_pred[:, 1, :]+1e-10) )
    #                     + SEQ_WEIGHT * y_true[:, 2, :] * torch.mul( torch.pow( torch.sub(1, y_pred[:, 2, :]), gamma ), torch.log(y_pred[:, 2, :]+1e-10) ))

    # This is cross entropy loss
    return - torch.mean(y_true[:, 0, :]*torch.log(y_pred[:, 0, :]+1e-10)
                        + y_true[:, 1, :]*torch.log(y_pred[:, 1, :]+1e-10)
                        + y_true[:, 2, :]*torch.log(y_pred[:, 2, :]+1e-10))


# def create_datapoints(seq, strand):

#     # seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
#     seq = seq.upper().replace('A', '1').replace('C', '2')
#     seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')
#     jn_start = JUNC_START
#     jn_end = JUNC_END
#     # print("Donor: ", seq[CL_MAX//2+jn_start: CL_MAX//2+jn_start+2])
#     # print("Donor: ", seq[CL_MAX//2+jn_end-2: CL_MAX//2+jn_end])

#     X0 = np.asarray(list(map(int, list(seq))))
#     Y0 = [np.zeros(SEQ_LEN) for t in range(1)]


#     if strand == '+':
#         for t in range(1):        
#             Y0[t][jn_start] = 2
#             Y0[t][jn_end] = 1
#     X, Y = one_hot_encode(X0, Y0)
#     return X, Y

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


def get_donor_acceptor_scores(D_YL, A_YL, D_YP, A_YP):
    # print("D_YL[:, 200]: ", D_YL[:, 200], "  ", len(D_YL[:, 200]))
    # print("D_YP[:, 200]: ", D_YP[:, 200], "  ", len(D_YP[:, 200]))

    return D_YL[:, 200], D_YP[:, 200], A_YL[:, 600], A_YP[:, 600]

def get_junc_scores(D_YL, A_YL, D_YP, A_YP, choice):
    # print("D_YL: ", D_YP[:, 200])
    # print("A_YL: ", A_YP[:, 600])
    if choice == "min":
        junc_labels = np.minimum(D_YL[:, 200], A_YL[:, 600])
        junc_scores = np.minimum(D_YP[:, 200], A_YP[:, 600])
    elif choice == "avg":
        junc_labels = np.minimum(D_YL[:, 200], A_YL[:, 600])
        junc_scores = np.mean([D_YP[:, 200], A_YP[:, 600]], axis=0)

    # print("junc_labels: ", junc_labels)
    # print("junc_scores: ", junc_scores)
    return junc_labels, junc_scores


def print_junc_statistics(D_YL, A_YL, D_YP, A_YP, threshold, TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN):

    label_junc_idx = (D_YL[:, 200]==1) & (A_YL[:, 600]==1)
    label_nonjunc_idx = (D_YL[:, 200]==0) & (A_YL[:, 600]==0)
    predict_junc_idx = (D_YP[:, 200]>=threshold) & (A_YP[:, 600]>=threshold)
    predict_nonjunc_idx = (D_YP[:, 200]<threshold) | (A_YP[:, 600]<threshold)

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



def junc_statistics(YL, YP, threshold, TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN):

    labels_1 = np.where(YL == 1)

    ####################
    # Donor 
    ####################
    thre_d = np.where(YP >= threshold)
    # print("thre_d: ", thre_d)
    # print("thre_a: ", thre_a)

    LCL_TOTAL_TP = len(np.intersect1d(labels_1, thre_d))
    LCL_TOTAL_FN = len(np.setdiff1d(labels_1, thre_d))
    LCL_TOTAL_FP = len(np.setdiff1d(thre_d, labels_1))
    LCL_TOTAL_TN = len(YL) - LCL_TOTAL_TP - LCL_TOTAL_FN - LCL_TOTAL_FP

    # print("Donor TPs: ", TPs)
    # print("Donor FNs: ", FNs)
    # print("Donor FPs: ", FPs)
    # print("Donor TNs: ", TNs)

    # LCL_TOTAL_TP = np.size(np.intersect1d(idx_true, idx_pred))
    # LCL_TOTAL_FN = len(idx_true) - LCL_TOTAL_TP
    # LCL_TOTAL_FP = len(idx_pred) - LCL_TOTAL_TP

    # # LCL_TOTAL_TN = np.size(np.intersect1d(label_nonjunc_idx, predict_nonjunc_idx))
    # LCL_TOTAL_TN = len(YL) - LCL_TOTAL_TP - LCL_TOTAL_FN - LCL_TOTAL_FP

    # print("LCL_TOTAL_TP: ", LCL_TOTAL_TP)
    # print("LCL_TOTAL_FN: ", LCL_TOTAL_FN)
    # print("LCL_TOTAL_FP: ", LCL_TOTAL_FP)
    # print("LCL_TOTAL_TN: ", LCL_TOTAL_TN)

    TOTAL_TP += LCL_TOTAL_TP
    TOTAL_FN += LCL_TOTAL_FN
    TOTAL_FP += LCL_TOTAL_FP
    TOTAL_TN += LCL_TOTAL_TN

    # # precision = np.size(np.intersect1d(idx_true, idx_pred)) \
    # #             / float(len(idx_pred))        
    # # recall = np.size(np.intersect1d(idx_true, idx_pred)) \
    # #             / float(len(idx_true)) 

    # # print("precision: ", precision)
    # # print("recall:    ", recall)
    return TOTAL_TP, TOTAL_FN, TOTAL_FP, TOTAL_TN, LCL_TOTAL_TP, LCL_TOTAL_FN, LCL_TOTAL_FP, LCL_TOTAL_TN