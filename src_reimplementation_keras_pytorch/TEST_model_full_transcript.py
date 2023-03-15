###############################################################################
# This file contains the code to train the SpliceAI model.
###############################################################################

from Informer2020.models.model import Informer
from SpliceNN_Informer import Informer_SpliceAI
import numpy as np
import sys
import time
import h5py
# import tensorflow.keras.backend as kb
import tensorflow as tf
from spliceai import *
# from SpliceNN import *
from spliceai_pytorch import *
from SpliceNN_dataset import *
# from SpliceNN_Informer import *
from utils_SpliceNN import *
from multi_gpu import *
from constants import *
from tqdm import tqdm
import math
import sys

###############################################################################
'''This parser takes as input the .h5 file produced by create_datafile.py and
outputs a .h5 file with datapoints of the form (X, Y), which can be understood
by Keras models.'''
###############################################################################


def model_fn(DNAs, labels, model, criterion, device):
    """Forward a batch through the model."""
    outs = model(DNAs)

    # print("labels: ", labels.size())
    # print("DNAs: ", DNAs.size())

    # labels = labels.sum(axis=1)
    # print("labels: ", labels.size())
    # print("outs: ", outs.size())

    loss = categorical_crossentropy_2d(labels, outs)
    return loss, outs

def categorical_crossentropy_2d(y_true, y_pred):
    # prod = output[:,0]*target
    # return -prod[prod<0].sum()
    # print("y_true: ", y_true)
    # print("Loss: ", - torch.mean(y_true[:, 0, :]*torch.log(y_pred[:, 0, :]+1e-10)
    #                     + y_true[:, 1, :]*torch.log(y_pred[:, 1, :]+1e-10)
    #                     + y_true[:, 2, :]*torch.log(y_pred[:, 2, :]+1e-10)))
    # print("y_true[:, 0, :]: ", y_true[:, 0, :])
    # print("y_pred[:, 0, :]: ", y_pred[:, 0, :])
    # print("y_true[:, 1, :]: ", y_true[:, 1, :])
    # print("y_pred[:, 1, :]: ", y_pred[:, 1, :])
    # print("y_true[:, 2, :]: ", y_true[:, 2, :])
    # print("y_pred[:, 2, :]: ", y_pred[:, 2, :])
    return - torch.mean(y_true[:, 0, :]*torch.log(y_pred[:, 0, :]+1e-10)
                        + y_true[:, 1, :]*torch.log(y_pred[:, 1, :]+1e-10)
                        + y_true[:, 2, :]*torch.log(y_pred[:, 2, :]+1e-10))

# fix random seed
def same_seeds(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

start_time = time.time()

assert sys.argv[1] in ['train', 'test', 'all']
assert sys.argv[2] in ['0', '1', 'all']

h5f = h5py.File(data_dir + 'datafile'
                + '_' + sys.argv[1] + '_' + sys.argv[2]
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

# h5f2 = h5py.File(data_dir + 'dataset'
#                 + '_' + sys.argv[1] + '_' + sys.argv[2]
#                 + '.h5', 'w')
BATCH_SIZE = 1
CHUNK_SIZE = 100

model_file = "MODEL/SpliceNN_spliceAI_pytroch_e_1059.pt"

start_time = time.time()

same_seeds(0)
device = torch.device("cuda" if torch.cuda.is_available() else "mps")
print(f"[Info]: Use {device} now!")

model = torch.load(model_file)
print(f"[Info]: Finish creating model!",flush = True)
print("model: ", model)

criterion = nn.CrossEntropyLoss()

start_time = time.time()



Y_true_1 = [[] for t in range(1)]
Y_true_2 = [[] for t in range(1)]
Y_pred_1 = [[] for t in range(1)]
Y_pred_2 = [[] for t in range(1)]

A_TOTAL_TP = [0]*8
A_TOTAL_FN = [0]*8
A_TOTAL_FP = [0]*8
A_TOTAL_TN = [0]*8

D_TOTAL_TP = [0]*8
D_TOTAL_FN = [0]*8
D_TOTAL_FP = [0]*8
D_TOTAL_TN = [0]*8
thresholds = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
TOTAL_TRANS = 0
TOTAL_A_SITES = 0
TOTAL_D_SITES = 0

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
        pbar = tqdm(total=len(X), ncols=0, desc="Train", unit=" step")

        test_dataset = myDataset(X, Y)
        test_loader = DataLoader(
            test_dataset,
            batch_size=1,
            shuffle=True,
            drop_last=True,
        )
        # print(">> In epoch: ", epoch_num, "  (current idx_train: ", i, ",  size: ", len(train_loader), ")")

        Yl = []
        Yp = []
        for batch_idx, data in enumerate(test_loader):
            # training_step += 1
            DNA, label = data 
            TOTAL_TRANS += label.shape[0]
            DNA = DNA.to(torch.float32).to(device)
            label = label.to(torch.float32).to(device)
            DNA = torch.permute(DNA, (0, 2, 1))
            label = torch.permute(label, (0, 2, 1))
            # print("\nDNA: ", DNA.size())
            # print("label: ", label.size())
            loss, yp = model_fn(DNA, label, model, criterion, device)
            is_expr = (label.sum(axis=(1,2)) >= 1)
            Acceptor_YL = label[is_expr, 1, :].flatten().to('cpu').detach().numpy()
            Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
            Donor_YL = label[is_expr, 2, :].flatten().to('cpu').detach().numpy()
            Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

            Yl.append(label)
            Yp.append(yp)
            for idx, threshold in enumerate(thresholds):
                print("threshold: ", threshold)
                A_TOTAL_TP[idx], A_TOTAL_FN[idx], A_TOTAL_FP[idx] = print_threshold_statistics(Acceptor_YL, Acceptor_YP, threshold, A_TOTAL_TP[idx], A_TOTAL_FN[idx], A_TOTAL_FP[idx])
                D_TOTAL_TP[idx], D_TOTAL_FN[idx], D_TOTAL_FP[idx] = print_threshold_statistics(Donor_YL, Donor_YP, threshold, D_TOTAL_TP[idx], D_TOTAL_FN[idx], D_TOTAL_FP[idx])

            # A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
            # D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
            # print(">> Acceptor: ")
            # print_topl_statistics(Acceptor_YL, Acceptor_YP)
            # print(">> Donor: ")
            # print_topl_statistics(Donor_YL, Donor_YP)

            batch_loss = loss.item()

            pbar.update(BATCH_SIZE)
            pbar.set_postfix(
                idx_train=len(test_loader),
                loss=f"{batch_loss:.6f}",
                # A_auc = f"{A_auc:.6f}",
                # D_auc = f"{D_auc:.6f}",
                # A_accuracy=f"{A_accuracy:.6f}",
                # D_accuracy=f"{D_accuracy:.6f}",
            )
            for idx, threshold in enumerate(thresholds):
                print(">> Threshold (", threshold, ")")
                print("\tAcceptor Precision: ", A_TOTAL_TP[idx]/(A_TOTAL_TP[idx]+A_TOTAL_FP[idx]))
                print("\tAcceptor Recall   : ", A_TOTAL_TP[idx]/(A_TOTAL_TP[idx]+A_TOTAL_FN[idx]))
                print("\tDonor Precision   : ", D_TOTAL_TP[idx]/(D_TOTAL_TP[idx]+D_TOTAL_FP[idx]))
                print("\tDonor Recall      : ", D_TOTAL_TP[idx]/(D_TOTAL_TP[idx]+D_TOTAL_FN[idx]))
                print("\tTotal testing transcripts   : ", TOTAL_TRANS)
                print("\tTotal labelled Acceptor splice sites : ", A_TOTAL_TP[idx] + A_TOTAL_FN[idx])
                print("\tTotal predicted Acceptor splice sites: ", A_TOTAL_TP[idx] + A_TOTAL_FP[idx])
                print("\tTotal labelled Donor splice sites    : ", D_TOTAL_TP[idx] + D_TOTAL_FN[idx])
                print("\tTotal predicted Donor splice sites   : ", D_TOTAL_TP[idx] + D_TOTAL_FP[idx])
                print("")
            # batch_accuracy = accuracy.item()
            # Updata model
        Yl = torch.stack(Yl, dim=1)
        Yp = torch.stack(Yp, dim=1)

        # for t in range(1):
        #     # print("Yl[t].sum(axis=(1,2)): ", len(Yl[t].sum(axis=(1,2))))
        #     is_expr = (Yl[t].sum(axis=(1,2)) >= 1)
        #     # print("is_expr: ", len(is_expr))
        #     Y_true_1[t].extend(Yl[t][is_expr, 1, :].flatten().to('cpu').detach().numpy())
        #     Y_true_2[t].extend(Yl[t][is_expr, 2, :].flatten().to('cpu').detach().numpy())
        #     Y_pred_1[t].extend(Yp[t][is_expr, 1, :].flatten().to('cpu').detach().numpy())
        #     Y_pred_2[t].extend(Yp[t][is_expr, 2, :].flatten().to('cpu').detach().numpy())

        # print("\n\n##########################################################")
        # print(">> valid_loader length: ", len(test_loader))
        # print("\n\033[1mValidation set metrics:\033[0m")
        # print("\n\033[1mAcceptor:\033[0m")
        # for t in range(1):
        #     print_topl_statistics(np.asarray(Y_true_1[t]),
        #                         np.asarray(Y_pred_1[t]))

        # print("\n\033[1mDonor:\033[0m")
        # for t in range(1):
        #     print_topl_statistics(np.asarray(Y_true_2[t]),
        #                         np.asarray(Y_pred_2[t]))
        # print("##########################################################\n\n\n")

    pbar.close()















    # X_batch = np.asarray(X_batch).astype('int8')
    # for t in range(1):
    #     Y_batch[t] = np.asarray(Y_batch[t]).astype('int8')


#     h5f2.create_dataset('X' + str(i), data=X_batch)
#     h5f2.create_dataset('Y' + str(i), data=Y_batch)

# h5f2.close()

###############################################################################
# Model
###############################################################################

# L = 32
# N_GPUS = 2

# if int(sys.argv[1]) == 80:
#     W = np.asarray([11, 11, 11, 11])
#     AR = np.asarray([1, 1, 1, 1])
#     BATCH_SIZE = 18*N_GPUS
# elif int(sys.argv[1]) == 400:
#     W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
#     AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
#     BATCH_SIZE = 18*N_GPUS
# elif int(sys.argv[1]) == 2000:
#     W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
#                     21, 21, 21, 21])
#     AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
#                      10, 10, 10, 10])
#     BATCH_SIZE = 12*N_GPUS
# elif int(sys.argv[1]) == 10000:
#     W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
#                     21, 21, 21, 21, 41, 41, 41, 41])
#     AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
#                      10, 10, 10, 10, 25, 25, 25, 25])

#     BATCH_SIZE = 6*N_GPUS

# Hyper-parameters:
# L: Number of convolution kernels
# W: Convolution window size in each residual unit
# AR: Atrous rate in each residual unit

# BATCH_SIZE = 24
# CL = 2 * np.sum(AR*(W-1))
# assert CL <= CL_max and CL == int(sys.argv[1])
# print("\033[1mContext nucleotides: %d\033[0m" % (CL))
# print("\033[1mSequence length (output): %d\033[0m" % (SL))

# model = SpliceAI(L, W, AR)
# model.summary()
# model_m = make_parallel(model, N_GPUS)
# model_m.compile(loss=categorical_crossentropy_2d, optimizer='adam')

###############################################################################
# Training and validation
###############################################################################

# def test_epoch():
    # idx = np.random.choice(idx_all)

    # for idx in idx_all:

# if __name__ == "__main__":
#     test_epoch()
# for epoch_num in range(EPOCH_NUM):
#     train_one_epoch(epoch_num)
#     torch.save(model, './MODELS_Conformer/SpliceNN_spliceAI_e_'+str(epoch_num)+'.pt')

























        # if training_step > 50:
        #     break
    # Yl = torch.stack(Yl, dim=1)
    # Yp = torch.stack(Yp, dim=1)

    # for t in range(1):
    #     is_expr = (Yl[t].sum(axis=(1,2)) >= 1)
    #     # print("is_expr: ", is_expr)
    #     # print("Yl[t]: ", Yl[t].size())
    #     # print("Yl[t]: ", Yl[t])
    #     # print("Yl[t].sum(axis=(1,2)): ", Yl[t].sum(axis=(1,2)))
    #     # print("Yl[t].sum(axis=(1,2)).shape: ", Yl[t].sum(axis=(1,2)).shape)
    #     Y_true_1[t].extend(Yl[t][is_expr, :, 1].flatten().to('cpu').detach().numpy())
    #     Y_true_2[t].extend(Yl[t][is_expr, :, 2].flatten().to('cpu').detach().numpy())
    #     Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten().to('cpu').detach().numpy())
    #     Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten().to('cpu').detach().numpy())
    # # CC += 1
    # # if CC == 3:
    # #     break
    # print("\n\n----------------------------------------------------------")
    # print(">> train_loader length: ", len(train_loader))
    # print("\033[1mTraining set metrics:\033[0m")
    # print("\n\033[1mAcceptor:\033[0m")
    # for t in range(1):
    #     print_topl_statistics(np.asarray(Y_true_1[t]),
    #                             np.asarray(Y_pred_1[t]))

    # print("\n\033[1mDonor:\033[0m")
    # for t in range(1):
    #     print_topl_statistics(np.asarray(Y_true_2[t]),
    #                             np.asarray(Y_pred_2[t]))
    # print("----------------------------------------------------------\n\n\n")







#     if (epoch_num+1) % len(idx_train) == 0:
#         # Printing metrics (see utils.py for details)
#         print("--------------------------------------------------------------")
#         print("\n\033[1mValidation set metrics:\033[0m")
#         Y_true_1 = [[] for t in range(1)]
#         Y_true_2 = [[] for t in range(1)]
#         Y_pred_1 = [[] for t in range(1)]
#         Y_pred_2 = [[] for t in range(1)]
#         for idx in idx_valid:
#             X = h5f['X' + str(idx)][:]
#             Y = h5f['Y' + str(idx)][:]
#             Xc, Yc = clip_datapoints(X, Y, 5000, N_GPUS) 
#             Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

#             if not isinstance(Yp, list):
#                 Yp = [Yp]

#             for t in range(1):

#                 is_expr = (Yc[t].sum(axis=(1,2)) >= 1)

#                 Y_true_1[t].extend(Yc[t][is_expr, :, 1].flatten())
#                 Y_true_2[t].extend(Yc[t][is_expr, :, 2].flatten())
#                 Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten())
#                 Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten())

#         print("\n\033[1mAcceptor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_1[t]),
#                                   np.asarray(Y_pred_1[t]))

#         print("\n\033[1mDonor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_2[t]),
#                                   np.asarray(Y_pred_2[t]))

#         print("\n\033[1mTraining set metrics:\033[0m")

#         Y_true_1 = [[] for t in range(1)]
#         Y_true_2 = [[] for t in range(1)]
#         Y_pred_1 = [[] for t in range(1)]
#         Y_pred_2 = [[] for t in range(1)]

#         for idx in idx_train[:len(idx_valid)]:

#             X = h5f['X' + str(idx)][:]
#             Y = h5f['Y' + str(idx)][:]

#             Xc, Yc = clip_datapoints(X, Y, CL, N_GPUS)
#             Yp = model_m.predict(Xc, batch_size=BATCH_SIZE)

#             if not isinstance(Yp, list):
#                 Yp = [Yp]

#             for t in range(1):

#                 is_expr = (Yc[t].sum(axis=(1,2)) >= 1)

#                 Y_true_1[t].extend(Yc[t][is_expr, :, 1].flatten())
#                 Y_true_2[t].extend(Yc[t][is_expr, :, 2].flatten())
#                 Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten())
#                 Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten())

#         print("\n\033[1mAcceptor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_1[t]),
#                                   np.asarray(Y_pred_1[t]))

#         print("\n\033[1mDonor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_2[t]),
#                                   np.asarray(Y_pred_2[t]))

#         print("Learning rate: %.5f" % (kb.get_value(model_m.optimizer.lr)))
#         print("--- %s seconds ---" % (time.time() - start_time))
#         start_time = time.time()

#         print("--------------------------------------------------------------")

#         model.save('./Models/SpliceAI' + sys.argv[1]
#                    + '_c' + sys.argv[2] + '.h5')

#         if (epoch_num+1) >= 6*len(idx_train):
#             kb.set_value(model_m.optimizer.lr,
#                          0.5*kb.get_value(model_m.optimizer.lr))
#             # Learning rate decay

# h5f.close()

































# TOTAL_UPDATE = 0
# for epoch_num in range(EPOCH_NUM):
#     for i, idx in enumerate(idx_train):
#         Y = h5f['Y' + str(idx)][:]
#         # print("Y.size(): ", Y.shape)
#         TOTAL_UPDATE += Y.shape[1]
# print(">> TOTAL_UPDATE: ", TOTAL_UPDATE)

# pbar = tqdm(total=TOTAL_UPDATE, ncols=0, desc="Train", unit=" step")
# training_step = 0
# for epoch_num in range(EPOCH_NUM):
#     ########################################
#     ## Training
#     ########################################
#     for i, idx in enumerate(idx_train):
#         X = h5f['X' + str(idx)][:]
#         Y = h5f['Y' + str(idx)][:]
#         Xc, Yc = clip_datapoints(X, Y, 5000, N_GPUS) 
#         # Xc = np.array(Xc)
#         # Yc = np.array(Yc)
#         # print("Xc.shape: ", Xc.shape)
#         # print("Yc.shape: ", Yc.shape)
#         train_dataset = myDataset(Xc, Yc)
#         train_loader = DataLoader(
#             train_dataset,
#             batch_size=BATCH_SIZE,
#             shuffle=True,
#             drop_last=True,
#             # num_workers=1,
#             # pin_memory=True
#             # collate_fn=collate_batch,
#         )
#         print(">> In epoch: ", epoch_num, "  (current idx_train: ", i, ",  size: ", len(train_loader), ")")

#         Y_true_1 = [[] for t in range(1)]
#         Y_true_2 = [[] for t in range(1)]
#         Y_pred_1 = [[] for t in range(1)]
#         Y_pred_2 = [[] for t in range(1)]
#         Yl = []
#         Yp = []
#         for data in train_loader:
#             training_step += 1
#             DNA, label = data 
#             # print("\nDNAs: ", DNAs.size())
#             # print("labels: ", labels.size())
#             DNA = DNA.to(torch.float32).to(device)
#             label = label.to(torch.float32).to(device)
#             loss, yp = model_fn(DNA, label, model, criterion, device)
#             # print("yp: ", yp)
#             Yl.append(label)
#             Yp.append(yp)
            
#             batch_loss = loss.item()

#             pbar.update()
#             pbar.set_postfix(
#                 loss=f"{batch_loss:.2f}",
#                 # accuracy=f"{batch_accuracy:.2f}",
#                 step = training_step,
#             )
#             # batch_accuracy = accuracy.item()
#             # Updata model
#             loss.backward()
#             optimizer.step()
#             scheduler.step()
#             optimizer.zero_grad()

#             if training_step > 50:
#                 break
#         Yl = torch.stack(Yl, dim=1)
#         Yp = torch.stack(Yp, dim=1)

#         for t in range(1):
#             is_expr = (Yl[t].sum(axis=(1,2)) >= 1)
#             # print("is_expr: ", is_expr)
#             # print("Yl[t]: ", Yl[t].size())
#             # print("Yl[t]: ", Yl[t])
#             # print("Yl[t].sum(axis=(1,2)): ", Yl[t].sum(axis=(1,2)))
#             # print("Yl[t].sum(axis=(1,2)).shape: ", Yl[t].sum(axis=(1,2)).shape)
#             Y_true_1[t].extend(Yl[t][is_expr, :, 1].flatten().to('cpu').detach().numpy())
#             Y_true_2[t].extend(Yl[t][is_expr, :, 2].flatten().to('cpu').detach().numpy())
#             Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten().to('cpu').detach().numpy())
#             Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten().to('cpu').detach().numpy())
#         # CC += 1
#         # if CC == 3:
#         #     break
#         print("\n\n----------------------------------------------------------")
#         print(">> train_loader length: ", len(train_loader))
#         print("\033[1mTraining set metrics:\033[0m")
#         print("\n\033[1mAcceptor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_1[t]),
#                                     np.asarray(Y_pred_1[t]))

#         print("\n\033[1mDonor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_2[t]),
#                                     np.asarray(Y_pred_2[t]))
#         print("----------------------------------------------------------\n\n\n")

#         if training_step > 20:
#             break

#     # ########################################
#     # ## Validation.
#     # ########################################
#     # TOTAL_VAL_UPDATE = 50
#     # counter_val = 0
#     # pbar_val = tqdm(total=TOTAL_VAL_UPDATE, ncols=0, desc="Train", unit=" step")

#     # Y_true_1 = [[] for t in range(1)]
#     # Y_true_2 = [[] for t in range(1)]
#     # Y_pred_1 = [[] for t in range(1)]
#     # Y_pred_2 = [[] for t in range(1)]

#     # for idx in idx_valid:
#     #     X = h5f['X' + str(idx)][:]
#     #     Y = h5f['Y' + str(idx)][:]
#     #     Xc, Yc = clip_datapoints(X, Y, 5000, N_GPUS) 
#     #     valid_dataset = myDataset(Xc, Yc)
#     #     valid_loader = DataLoader(
#     #         valid_dataset,
#     #         batch_size=BATCH_SIZE,
#     #         shuffle=True,
#     #     )
#     #     Yl = []
#     #     Yp = []
#     #     for data in valid_loader:
#     #         counter_val += 1
#     #         DNA, label = data 
#     #         # print("\nDNAs: ", DNAs.size())
#     #         # print("labels: ", labels.size())
#     #         DNA = DNA.to(torch.float32).to(device)
#     #         label = label.to(torch.float32).to(device)
#     #         loss, yp = model_fn(DNA, label, model, criterion, device)
#     #         Yl.append(label)
#     #         Yp.append(yp)
#     #         batch_loss = loss.item()
#     #         pbar_val.update()
#     #         pbar_val.set_postfix(
#     #             loss=f"{batch_loss:.2f}",
#     #             step = training_step,
#     #         )
#     #         if counter_val >= TOTAL_VAL_UPDATE:
#     #             pbar_val.close()
#     #             break

#     #     Yl = torch.stack(Yl, dim=1)
#     #     Yp = torch.stack(Yp, dim=1)

#     #     for t in range(1):
#     #         # print("Yl[t].sum(axis=(1,2)): ", len(Yl[t].sum(axis=(1,2))))
#     #         is_expr = (Yl[t].sum(axis=(1,2)) >= 1)
#     #         # print("is_expr: ", len(is_expr))
#     #         Y_true_1[t].extend(Yl[t][is_expr, :, 1].flatten().to('cpu').detach().numpy())
#     #         Y_true_2[t].extend(Yl[t][is_expr, :, 2].flatten().to('cpu').detach().numpy())
#     #         Y_pred_1[t].extend(Yp[t][is_expr, :, 1].flatten().to('cpu').detach().numpy())
#     #         Y_pred_2[t].extend(Yp[t][is_expr, :, 2].flatten().to('cpu').detach().numpy())

#     #     print("\n\n##########################################################")
#     #     print(">> valid_loader length: ", len(valid_loader))
#     #     print("\n\033[1mValidation set metrics:\033[0m")
#     #     print("\n\033[1mAcceptor:\033[0m")
#     #     for t in range(1):
#     #         print_topl_statistics(np.asarray(Y_true_1[t]),
#     #                             np.asarray(Y_pred_1[t]))

#     #     print("\n\033[1mDonor:\033[0m")
#     #     for t in range(1):
#     #         print_topl_statistics(np.asarray(Y_true_2[t]),
#     #                             np.asarray(Y_pred_2[t]))
#     #     print("##########################################################\n\n\n")
#     #     if counter_val >= TOTAL_VAL_UPDATE:
#     #         break

#     torch.save(model, './MODEL/SpliceNN_e_'+str(epoch_num)+'.pt')

#         # if (epoch_num+1) >= 6*len(idx_train):
#         #     kb.set_value(model_m.optimizer.lr,
#         #                  0.5*kb.get_value(model_m.optimizer.lr))





#             # Learning rate decay

#     # try:
#     #     batch = next(train_iterator)
#     # except StopIteration:
#     #     train_iterator = iter(train_loader)
#     #     batch = next(train_iterator)
    
#     # torch.FloatTensor(speaker).long()




h5f.close()
        
###############################################################################

