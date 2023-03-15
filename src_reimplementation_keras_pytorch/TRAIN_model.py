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
from SPLAM import *
from spliceai_pytorch import *
from SpliceNN_dataset import *
# from SpliceNN_Informer import *
from utils_SpliceNN import *
from multi_gpu import *
from constants import *

# assert int(sys.argv[1]) in [80, 400, 2000, 10000]

###############################################################################
# Model
###############################################################################

L = 32
N_GPUS = 2

BATCH_SIZE = 2
N_WORKERS = 8
# Hyper-parameters:
# L: Number of convolution kernels
# W: Convolution window size in each residual unit
# AR: Atrous rate in each residual unit
# print("\033[1mContext nucleotides: %d\033[0m" % (CL))
print("\033[1mSequence length (output): %d\033[0m" % (SL))

# model = SpliceAI(L, W, AR)
# model.summary()
# model_m = make_parallel(model, N_GPUS)
# model_m.compile(loss=categorical_crossentropy_2d, optimizer='adam')

###############################################################################
# Training and validation
###############################################################################

h5f = h5py.File(data_dir + 'data/dataset_train_all.h5', 'r')
print(h5f.keys())
num_idx = len(list(h5f.keys()))//2
idx_all = np.random.permutation(num_idx)
idx_train = idx_all[:int(0.9*num_idx)]
idx_valid = idx_all[int(0.9*num_idx):]

print("len(idx_train): ", len(idx_train))
print("idx_train: ", idx_train)
print("idx_valid: ", idx_valid)

np.random.shuffle(idx_train)
np.random.shuffle(idx_valid)

EPOCH_NUM = 10
start_time = time.time()

i = "config1"
same_seeds(0)
device = torch.device("cuda" if torch.cuda.is_available() else "mps")
print(f"[Info]: Use {device} now!")


model_config = train_parse_args()["model_config"][i]
model = SpliceNN(model_config, seq_len=5000).to(device)
criterion = nn.CrossEntropyLoss()
optimizer = AdamW(model.parameters(), lr=1e-3)
scheduler = get_cosine_schedule_with_warmup(optimizer, train_parse_args()["warmup_steps"], train_parse_args()["total_steps"])
print(f"[Info]: Finish creating model!",flush = True)
print("model: ", model)



EPOCH_NUM = 10*len(idx_train)
start_time = time.time()


def train_one_epoch(epoch_idx):

    idx = np.random.choice(idx_train)

    X = h5f['X' + str(idx)][:]
    Y = h5f['Y' + str(idx)][:]

    Xc, Yc = clip_datapoints(X, Y, 5000, N_GPUS) 
    pbar = tqdm(total=len(Xc), ncols=0, desc="Train", unit=" step")

    train_dataset = myDataset(Xc, Yc)
    train_loader = DataLoader(
        train_dataset,
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=True,
        # num_workers=1,
        # pin_memory=True
        # collate_fn=collate_batch,
    )
    print(">> In epoch: ", epoch_num, "  (current idx_train: ", i, ",  size: ", len(train_loader), ")")

    for batch_idx, data in enumerate(train_loader):
        # training_step += 1
        DNA, label = data 
        # print("\nDNAs: ", DNAs.size())
        # print("labels: ", labels.size())
        DNA = DNA.to(torch.float32).to(device)
        label = label.to(torch.float32).to(device)
        loss, yp = model_fn(DNA, label, model, criterion, device)
        is_expr = (label.sum(axis=(1,2)) >= 1)
        Acceptor_YL = label[is_expr, :, 1].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, :, 1].flatten().to('cpu').detach().numpy()
        Donor_YL = label[is_expr, :, 2].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, :, 2].flatten().to('cpu').detach().numpy()

        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        batch_loss = loss.item()

        pbar.update(BATCH_SIZE)
        pbar.set_postfix(
            epoch=epoch_num,
            idx_train=len(train_loader),
            loss=f"{batch_loss:.6f}",
            A_accuracy=f"{A_accuracy:.6f}",
            D_accuracy=f"{D_accuracy:.6f}",
            A_auc = f"{A_auc:.6f}",
            D_auc = f"{D_auc:.6f}",
        )
        # batch_accuracy = accuracy.item()
        # Updata model
        loss.backward()
        optimizer.step()
        scheduler.step()
        optimizer.zero_grad()    


for epoch_num in range(EPOCH_NUM):
    train_one_epoch(epoch_num)
    torch.save(model, './MODELS_Conformer/SpliceNN_spliceAI_e_'+str(epoch_num)+'.pt')


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




pbar.close()
h5f.close()
        
###############################################################################

