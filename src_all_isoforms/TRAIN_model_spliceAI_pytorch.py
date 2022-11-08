###############################################################################
# This file contains the code to train the SpliceAI model.
###############################################################################

import numpy as np
import sys
import time
import h5py
# import tensorflow.keras.backend as kb
import tensorflow as tf
# from spliceai import *
# from SpliceNN import *
# from SpliceNN_dataset import *
from spliceai_pytorch import *
# from SpliceNN_Informer import *
from SpliceNN_dataset import *
from utils import *
from multi_gpu import *
from constants import *
from utils_SpliceNN import *  
from torch.optim import Optimizer, AdamW
from torch.optim.lr_scheduler import LambdaLR
from tqdm import tqdm
import math

L = 32
N_GPUS = 2
CL = 0
BATCH_SIZE = 40
EPOCH_NUM = 10
h5f = h5py.File(data_dir + 'data/dataset_train.h5', 'r')
device = torch.device("cuda" if torch.cuda.is_available() else "mps")
print(f"[Info]: Use {device} now!")

assert int(sys.argv[1]) in [80, 400, 2000, 10000]

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

def model_fn(DNAs, labels, model, criterion, device):
    """Forward a batch through the model."""
    outs = model(DNAs)

    # print("labels: ", labels.size())
    # print("DNAs: ", DNAs.size())

    # labels = labels.sum(axis=1)
    # print("labels: ", labels.size())

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

###############################################################################
# Model
###############################################################################

def train_one_epoch(epoch_idx, idx_train, model, criterion, optimizer):
    idx = np.random.choice(idx_train)

    X = h5f['X' + str(idx)][:]
    Y = h5f['Y' + str(idx)][:]

    print("CL: ", CL)
    Xc, Yc = clip_datapoints_spliceAI(X, Y, CL, N_GPUS)
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
    for batch_idx, data in enumerate(train_loader):
        # training_step += 1
        DNA, label = data 
        # print("\nDNAs: ", DNAs.size())
        # print("labels: ", labels.size())
        DNA = DNA.to(torch.float32).to(device)
        label = label.to(torch.float32).to(device)
        DNA = torch.permute(DNA, (0, 2, 1))
        label = torch.permute(label, (0, 2, 1))
        loss, yp = model_fn(DNA, label, model, criterion, device)
        is_expr = (label.sum(axis=(1,2)) >= 1)
        # print("is_expr: ", is_expr)

        Acceptor_YL = label[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = label[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        batch_loss = loss.item()

        pbar.update(BATCH_SIZE)
        pbar.set_postfix(
            epoch=epoch_idx,
            idx_train=len(train_loader),
            loss=f"{batch_loss:.6f}",
            A_accuracy=f"{A_accuracy:.6f}",
            D_accuracy=f"{D_accuracy:.6f}",
            A_auc = f"{A_auc:.6f}",
            D_auc = f"{D_auc:.6f}",
        )
        loss.backward()
        optimizer.step()
        # scheduler.step()
        optimizer.zero_grad()
    # if (epoch_num+1) % len(idx_train) == 1:
    #     break
    # if True:
    pbar.close()
    
    
# def test_one_epoch(epoch_idx):
#     print("--------------------------------------------------------------")
#     print("\n\033[1mValidation set metrics:\033[0m")

#     Y_true_1 = [[] for t in range(1)]
#     Y_true_2 = [[] for t in range(1)]
#     Y_pred_1 = [[] for t in range(1)]
#     Y_pred_2 = [[] for t in range(1)]

#     for idx in idx_valid[1:4]:
#         X = h5f['X' + str(idx)][:]
#         Y = h5f['Y' + str(idx)][:]
#         Xc, Yc = clip_datapoints_spliceAI(X, Y, CL, N_GPUS)
#         pbar_val = tqdm(total=len(Xc), ncols=0, desc="Train", unit=" step")
#         valid_dataset = myDataset(Xc, Yc)
#         valid_loader = DataLoader(
#             valid_dataset,
#             batch_size=BATCH_SIZE,
#             shuffle=True,
#         )
#         Yl = []
#         Yp = []
#         for data in valid_loader:
#             # counter_val += 1
#             DNA, label = data 
#             # print("\nDNAs: ", DNAs.size())
#             # print("labels: ", labels.size())
#             DNA = DNA.to(torch.float32).to(device)
#             label = label.to(torch.float32).to(device)
#             # print("\nDNAs: ", DNA.size())
#             # print("labels: ", label.size())
#             DNA = torch.permute(DNA, (0, 2, 1))
#             label = torch.permute(label, (0, 2, 1))
#             loss, yp = model_fn(DNA, label, model, criterion, device)
#             # is_expr = (label.sum(axis=(1,2)) >= 1)


#             Yl.append(label)
#             Yp.append(yp)
#             batch_loss = loss.item()
#             pbar_val.update(BATCH_SIZE)
#             pbar_val.set_postfix(
#                 loss=f"{batch_loss:.2f}",
#             )
#             # if counter_val >= TOTAL_VAL_UPDATE:
#             #     pbar_val.close()
#             #     break

#         Yl = torch.stack(Yl, dim=1)
#         Yp = torch.stack(Yp, dim=1)

#         for t in range(1):
#             # print("Yl[t].sum(axis=(1,2)): ", len(Yl[t].sum(axis=(1,2))))
#             is_expr = (Yl[t].sum(axis=(1,2)) >= 1)
#             # print("is_expr: ", len(is_expr))
#             Y_true_1[t].extend(Yl[t][is_expr, 1, :].flatten().to('cpu').detach().numpy())
#             Y_true_2[t].extend(Yl[t][is_expr, 2, :].flatten().to('cpu').detach().numpy())
#             Y_pred_1[t].extend(Yp[t][is_expr, 1, :].flatten().to('cpu').detach().numpy())
#             Y_pred_2[t].extend(Yp[t][is_expr, 2, :].flatten().to('cpu').detach().numpy())

#         print("\n\n##########################################################")
#         print(">> valid_loader length: ", len(valid_loader))
#         print("\n\033[1mValidation set metrics:\033[0m")
#         print("\n\033[1mAcceptor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_1[t]),
#                                 np.asarray(Y_pred_1[t]))

#         print("\n\033[1mDonor:\033[0m")
#         for t in range(1):
#             print_topl_statistics(np.asarray(Y_true_2[t]),
#                                 np.asarray(Y_pred_2[t]))
#         print("##########################################################\n\n\n")


def main():
    global CL
    if int(sys.argv[1]) == 80:
        W = np.asarray([11, 11, 11, 11])
        AR = np.asarray([1, 1, 1, 1])
        BATCH_SIZE = 18*N_GPUS
    elif int(sys.argv[1]) == 400:
        W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
        AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
        BATCH_SIZE = 18*N_GPUS
    elif int(sys.argv[1]) == 2000:
        W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                        21, 21, 21, 21])
        AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                        10, 10, 10, 10])
        BATCH_SIZE = 12*N_GPUS
    elif int(sys.argv[1]) == 10000:
        W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                        21, 21, 21, 21, 41, 41, 41, 41])
        AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                        10, 10, 10, 10, 25, 25, 25, 25])

        BATCH_SIZE = 6*N_GPUS

    # Hyper-parameters:
    # L: Number of convolution kernels
    # W: Convolution window size in each residual unit
    # AR: Atrous rate in each residual unit


    CL = 2 * np.sum(AR*(W-1))
    assert CL <= CL_max and CL == int(sys.argv[1])
    print("\033[1mContext nucleotides: %d\033[0m" % (CL))
    print("\033[1mSequence length (output): %d\033[0m" % (SL))

    # h5f = h5py.File(data_dir + 'data/dataset_train.h5', 'r')

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

    start_time = time.time()

    i = "config1"
    same_seeds(0)


    model = SpliceAI(L, W, AR).to(device)
    criterion = nn.CrossEntropyLoss()
    optimizer = AdamW(model.parameters(), lr=1e-4)
    # scheduler = get_cosine_schedule_with_warmup(optimizer, 1000, 200000)
    print(f"[Info]: Finish creating model!",flush = True)
    print("model: ", model)


    EPOCH_NUM = 10*len(idx_train)
    start_time = time.time()

    for epoch_num in range(EPOCH_NUM):
        train_one_epoch(epoch_num, idx_train, model, criterion, optimizer)
        if epoch_num % 10 == 0:
            torch.save(model, './MODEL/SpliceNN_spliceAI_pytroch_e_'+str(epoch_num)+'.pt')

if __name__ == "__main__":
    main()

###############################################################################