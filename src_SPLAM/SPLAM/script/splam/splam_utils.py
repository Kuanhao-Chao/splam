###############################################################################
'''This code has functions which process the information in the .h5 files
datafile_{}_{}.h5 and convert them into a format usable by Keras.'''
###############################################################################

import torch
from torch.optim.lr_scheduler import LambdaLR
from torch.optim import Optimizer
from torch.utils.data import Dataset, DataLoader, random_split

import numpy as np
import math
import random
import pickle
from sklearn.metrics import average_precision_score

SEQ_LEN = 800
JUNC_START = 200
JUNC_END = 600
CL_MAX = 10000
SL = 800
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

           
def one_hot_encode_classifier(Xd):
    return IN_MAP[Xd.astype('int8')]
           
#######################################
# This is for Conformer model 
#######################################
def create_datapoints(seq, strand):
    # seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0').replace('K', '0').replace('R', '0').replace('Y', '0').replace('M', '0')
    jn_start = JUNC_START
    jn_end = JUNC_END

    #######################################
    # predicting splice / non-splice
    #######################################
    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = 0
    if strand == '+':
        Y0 = 1
    X = one_hot_encode_classifier(X0)
    return X, Y0

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

def get_accuracy(y_prob, y_true):
    assert y_true.ndim == 1 and y_true.size() == y_prob.size()
    y_prob = y_prob > 0.5
    return (y_true == y_prob).sum().item() / y_true.size(0)

def model_fn(DNAs, labels, model, criterion):
    """Forward a batch through the model."""
    outs = model(DNAs)
    outs = torch.flatten(outs)
    loss, accuracy = categorical_crossentropy_2d(labels, outs, criterion)
    return loss, accuracy, outs

def weighted_binary_cross_entropy(output, target, weights=None):    
    if weights is not None:
        assert len(weights) == 2
        loss = weights[1] * (target * torch.log(output+1e-10)) + \
               weights[0] * ((1 - target) * torch.log(1 - output+1e-10))
    else:
        loss = target * torch.log(output+1e-10) + (1 - target) * torch.log(1 - output+1e-10)
    return torch.neg(torch.mean(loss))

def categorical_crossentropy_2d(y_true, y_pred, criterion):
    # WEIGHT = 800
    # print("y_true: ", y_true)
    # print("y_pred: ", y_pred)
    weights = torch.FloatTensor([1.0, 1.0]) 
    return weighted_binary_cross_entropy(y_pred, y_true, weights), get_accuracy(y_pred, y_true)
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

    # return - torch.mean(y_true[:, 0, :]*torch.log(y_pred[:, 0, :]+1e-10)
    #                     + WEIGHT*y_true[:, 1, :]*torch.log(y_pred[:, 1, :]+1e-10)
    #                     + WEIGHT*y_true[:, 2, :]*torch.log(y_pred[:, 2, :]+1e-10))

def print_top_1_statistics(y_true, y_pred):
    
    idx_true = np.nonzero(y_true == 1)[0]
    argsorted_y_pred = np.argsort(y_pred)
    sorted_y_pred = np.sort(y_pred)

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






























# # Random_90_10 / Chromosome_90_10
# # TARGET = "Chromosome_split_p_n_nn_n1"
# SEQ_LEN = "800"
# # os.makedirs("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/", exist_ok=True)

# def split_seq_name(seq):
#     return seq[1:]



# class myDataset(Dataset):
#     def __init__(self, type, of, shuffle, segment_len=800):
#         self.segment_len = segment_len
#         self.data = []
#         self.indices = []
#         if type == "train":
#             CONSTANT_SIZE = 176086
#         else:
#             CONSTANT_SIZE = 23914

#         # CONSTANT_SIZE = 500
#         CONSTANT_SIZE_NEG = math.ceil(CONSTANT_SIZE*2/3)
#         #################################
#         ## Processing 'POSITIVE' samples
#         #################################
#         pidx = 0


#         with open(of, "r") as f:
#             print("of: ", of)
#             lines = f.read().splitlines()
#             seq_name = ""
#             seq = ""
#             for line in lines:
#                 # print(line)
#                 if pidx % 2 == 0:
#                     seq_name = split_seq_name(line)
#                 elif pidx % 2 == 1:
#                     seq = line
#                     if seq[0] == ">":
#                         seq_name = line
#                         continue
                    
#                     X, Y = create_datapoints(seq, '+')
#                     X = torch.Tensor(np.array(X))
#                     # print(X)
#                     # print(Y)
#                     if X.size()[0] != 800:
#                         print("seq_name: ", seq_name)
#                         print(X.size())
#                         # print(Y.size())
#                     self.data.append([X, Y, seq_name])
#                 pidx += 1
#                 if pidx %10000 == 0:
#                     print("pidx: ", pidx)

#         index_shuf = list(range(len(self.data)))
#         if shuffle:
#             random.shuffle(index_shuf)
#             # Shuffle just in a certain range.

#         list_shuf = [self.data[i] for i in index_shuf]
#         self.data = list_shuf 
#         self.indices = index_shuf
#         print("pidx: ", pidx)

#     def __len__(self):
#         return len(self.data)
 
#     def __getitem__(self, index):
#         # Load preprocessed mel-spectrogram.
#         # print("self.data: ", self.data[index])
#         feature = self.data[index][0]
#         label = self.data[index][1]
#         seq_name = self.data[index][2]

#         feature = torch.flatten(feature, start_dim=1)
#         return feature, label, seq_name



# def get_dataloader(batch_size, n_workers, output_file, shuffle):
#     """Generate dataloader"""
#     # testset = myDataset("test", output_file, int(SEQ_LEN))
#     testset = myDataset("test", output_file, shuffle, int(SEQ_LEN))

#     print("testset.indices: ", testset.indices)

#     test_loader = DataLoader(
#         testset,
#         batch_size = batch_size,
#         # batch_size=len(validset),
#         # num_workers=n_workers,
#         drop_last=True,
#         pin_memory=True,
#         # collate_fn=collate_batch,
#     )
#     if batch_size == 1:
#         print("shuffle: ", shuffle)
#         torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.nobatch.pt")
#         with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.nobatch.indices.pkl", "wb") as f:
#             pickle.dump(testset.indices, f)
            
#     elif shuffle:
#         print("shuffle: ", shuffle)
#         torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.shuffle.pt")
#         with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.shuffle.indices.pkl", "wb") as f:
#             pickle.dump(testset.indices, f)

#     else:
#         print("shuffle: ", shuffle)
#         torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.noshuffle.pt")
#         with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.noshuffle.indices.pkl", "wb") as f:
#             pickle.dump(testset.indices, f)

#     return test_loader

# # def get_dataloader(batch_size, n_workers):
# #     # print("Loading dataset: ", "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
# #     train_loader = torch.load("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")

# #     print("Loading dataset: ", "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
# #     test_loader = torch.load("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
# #     # return test_loader
# #     return train_loader, test_loader

# def get_dataloader(batch_size, n_workers, output_file, shuffle, repeat_idx):
#     """Generate dataloader"""
#     # testset = myDataset("test", output_file, int(SEQ_LEN))
#     print("output_file: ", output_file)
#     testset = myDataset("test", output_file, shuffle, int(SEQ_LEN))

#     # print("testset.indices: ", testset.indices)

#     test_loader = DataLoader(
#         testset,
#         batch_size = batch_size,
#         # batch_size=len(validset),
#         # num_workers=n_workers,
#         drop_last=True,
#         pin_memory=True,
#         # collate_fn=collate_batch,
#     )
#     if batch_size == 1:
#         print("shuffle: ", shuffle)
#         # torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.nobatch."+str(repeat_idx)+"."+str(batch_size)+".pt")
#         with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.nobatch.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
#             pickle.dump(testset.indices, f)
            
#     elif shuffle:
#         print("shuffle: ", shuffle)
#         # torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.shuffle."+str(repeat_idx)+"."+str(batch_size)+".pt")
#         with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.shuffle.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
#             pickle.dump(testset.indices, f)

#     else:
#         print("shuffle: ", shuffle)
#         # torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.noshuffle."+str(repeat_idx)+"."+str(batch_size)+".pt")
#         with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.noshuffle.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
#             pickle.dump(testset.indices, f)

#     return test_loader