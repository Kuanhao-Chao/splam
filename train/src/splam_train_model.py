import torch
import torch.nn as nn
from torch.optim import Optimizer, AdamW
from torch.optim.lr_scheduler import LambdaLR
import numpy as np
from tqdm import tqdm
import os
import platform
from sklearn.model_selection import KFold

from splam_utils import *
from splam_dataset_Chromsome import *
from SPLAM import *

#############################
# Global variable definition
#############################
EPOCH_NUM = 15
BATCH_SIZE = 100
N_WORKERS = 1
L = 64
SEQ_LEN = "800"
JUNC_THRESHOLD = 0.1

W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                11, 11, 11, 11, 21, 21, 21, 21,
                21, 21, 21, 21])
AR = np.asarray([1, 1, 1, 1, 5, 5, 5, 5,
                 10, 10, 10, 10, 15, 15, 15, 15,
                20, 20, 20, 20])

CL = 2 * np.sum(AR*(W-1))
print("\033[1mSequence length (input nucleotides): %d\033[0m" % (SL))
same_seeds(0)

#############################
# Selecting device
#############################
device_str = None
if torch.cuda.is_available():
    device_str = "cuda"
else:
    if platform.system() == "Darwin":
        device_str = "mps"
    else:
        device_str = "cpu"
device = torch.device(device_str)
print(f"\033[1m[Info]: Use {device} now!\033[0m")


#############################
# Model Initialization
#############################
model = SPLAM(L, W, AR).to(device)
criterion = nn.BCELoss()
optimizer = AdamW(model.parameters(), lr=1e-3)
print(f"\033[1m[Info]: Finish creating model!\033[0m",flush = True)
print("######################################################")
print("Model: ", model)
print("######################################################\n")

#############################
# Creating directories
#############################
MODEL_VERSION = "SPLAM_v13/"
MODEL_OUTPUT_BASE = "./MODEL/"+MODEL_VERSION
LOG_OUTPUT_BASE = MODEL_OUTPUT_BASE + "LOG/"
LOG_OUTPUT_TRAIN_BASE = MODEL_OUTPUT_BASE + "LOG/TRAIN/"
LOG_OUTPUT_VAL_BASE = MODEL_OUTPUT_BASE + "LOG/VAL/"
LOG_OUTPUT_TEST_BASE = MODEL_OUTPUT_BASE + "LOG/TEST/"
os.makedirs("./INPUTS/"+SEQ_LEN+"bp/"+MODEL_VERSION, exist_ok=True)
os.makedirs(LOG_OUTPUT_TRAIN_BASE, exist_ok=True)
os.makedirs(LOG_OUTPUT_VAL_BASE, exist_ok=True)
os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)

#############################
# Training Data initialization
#############################
# save_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS)
# train_loader, val_loader, test_loader = get_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS)
train_loader, test_loader = get_dataloader(BATCH_SIZE, MODEL_VERSION, N_WORKERS)
train_iterator = iter(train_loader)
# valid_iterator = iter(val_loader)
test_iterator = iter(test_loader)
print(f"\033[1m[Info]: Finish loading data!\033[0m",flush = True)
print("train_iterator: ", len(train_loader))
# print("valid_iterator: ", len(val_loader))
print("valid_iterator: ", len(test_loader))



#############################
# Initialize scheduler
#############################
scheduler = get_cosine_schedule_with_warmup(optimizer, 1000, len(train_loader)*EPOCH_NUM)
print(f"[Info]: Initialized the scheduler! Warmup steps: ", 1000, ";  Total steps: ", len(train_loader)*EPOCH_NUM)

############################
# Log for training
############################
train_log_loss = LOG_OUTPUT_TRAIN_BASE + "train_loss.txt"
train_log_acc = LOG_OUTPUT_TRAIN_BASE + "train_accuracy.txt"
train_log_lr = LOG_OUTPUT_TRAIN_BASE + "train_lr.txt"

train_log_A_auprc = LOG_OUTPUT_TRAIN_BASE + "train_A_auprc.txt"
train_log_A_threshold_precision = LOG_OUTPUT_TRAIN_BASE + "train_A_threshold_precision.txt"
train_log_A_threshold_recall = LOG_OUTPUT_TRAIN_BASE + "train_A_threshold_recall.txt"
train_log_D_auprc = LOG_OUTPUT_TRAIN_BASE + "train_D_auprc.txt"
train_log_D_threshold_precision = LOG_OUTPUT_TRAIN_BASE + "train_D_threshold_precision.txt"
train_log_D_threshold_recall = LOG_OUTPUT_TRAIN_BASE + "train_D_threshold_recall.txt"
train_log_J_threshold_precision = LOG_OUTPUT_TRAIN_BASE + "train_J_threshold_precision.txt"
train_log_J_threshold_recall = LOG_OUTPUT_TRAIN_BASE + "train_J_threshold_recall.txt"

fw_train_log_loss = open(train_log_loss, 'w')
fw_train_log_acc = open(train_log_acc, 'w')
fw_train_log_lr = open(train_log_lr, 'w')

fw_train_log_A_auprc = open(train_log_A_auprc, 'w')
fw_train_log_A_threshold_precision = open(train_log_A_threshold_precision, 'w')
fw_train_log_A_threshold_recall = open(train_log_A_threshold_recall, 'w')
fw_train_log_D_auprc = open(train_log_D_auprc, 'w')
fw_train_log_D_threshold_precision = open(train_log_D_threshold_precision, 'w')
fw_train_log_D_threshold_recall = open(train_log_D_threshold_recall, 'w')
fw_train_log_J_threshold_precision = open(train_log_J_threshold_precision, 'w')
fw_train_log_J_threshold_recall = open(train_log_J_threshold_recall, 'w')

############################
# Log for validation
############################
val_log_loss = LOG_OUTPUT_VAL_BASE + "val_loss.txt"
val_log_acc = LOG_OUTPUT_VAL_BASE + "val_accuracy.txt"

val_log_A_auprc = LOG_OUTPUT_VAL_BASE + "val_A_auprc.txt"
val_log_A_threshold_precision = LOG_OUTPUT_VAL_BASE + "val_A_threshold_precision.txt"
val_log_A_threshold_recall = LOG_OUTPUT_VAL_BASE + "val_A_threshold_recall.txt"
val_log_D_auprc = LOG_OUTPUT_VAL_BASE + "val_D_auprc.txt"
val_log_D_threshold_precision = LOG_OUTPUT_VAL_BASE + "val_D_threshold_precision.txt"
val_log_D_threshold_recall = LOG_OUTPUT_VAL_BASE + "val_D_threshold_recall.txt"
val_log_J_threshold_precision = LOG_OUTPUT_VAL_BASE + "val_J_threshold_precision.txt"
val_log_J_threshold_recall = LOG_OUTPUT_VAL_BASE + "val_J_threshold_recall.txt"

fw_val_log_loss = open(val_log_loss, 'w')
fw_val_log_acc = open(val_log_acc, 'w')

fw_val_log_A_auprc = open(val_log_A_auprc, 'w')
fw_val_log_A_threshold_precision = open(val_log_A_threshold_precision, 'w')
fw_val_log_A_threshold_recall = open(val_log_A_threshold_recall, 'w')
fw_val_log_D_auprc = open(val_log_D_auprc, 'w')
fw_val_log_D_threshold_precision = open(val_log_D_threshold_precision, 'w')
fw_val_log_D_threshold_recall = open(val_log_D_threshold_recall, 'w')
fw_val_log_J_threshold_precision = open(val_log_J_threshold_precision, 'w')
fw_val_log_J_threshold_recall = open(val_log_J_threshold_recall, 'w')

############################
# Log for testing
############################
test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
test_log_acc = LOG_OUTPUT_TEST_BASE + "test_accuracy.txt"

test_log_A_auprc = LOG_OUTPUT_TEST_BASE + "test_A_auprc.txt"
test_log_A_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_A_threshold_precision.txt"
test_log_A_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_A_threshold_recall.txt"
test_log_D_auprc = LOG_OUTPUT_TEST_BASE + "test_D_auprc.txt"
test_log_D_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_D_threshold_precision.txt"
test_log_D_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_D_threshold_recall.txt"
test_log_J_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_J_threshold_precision.txt"
test_log_J_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_J_threshold_recall.txt"

fw_test_log_loss = open(test_log_loss, 'w')
fw_test_log_acc = open(test_log_acc, 'w')

fw_test_log_A_auprc = open(test_log_A_auprc, 'w')
fw_test_log_A_threshold_precision = open(test_log_A_threshold_precision, 'w')
fw_test_log_A_threshold_recall = open(test_log_A_threshold_recall, 'w')
fw_test_log_D_auprc = open(test_log_D_auprc, 'w')
fw_test_log_D_threshold_precision = open(test_log_D_threshold_precision, 'w')
fw_test_log_D_threshold_recall = open(test_log_D_threshold_recall, 'w')
fw_test_log_J_threshold_precision = open(test_log_J_threshold_precision, 'w')
fw_test_log_J_threshold_recall = open(test_log_J_threshold_recall, 'w')


def get_lr(optimizer):
    for param_group in optimizer.param_groups:
        return param_group['lr']


def train_one_epoch(epoch_idx, train_loader):
    epoch_loss = 0
    epoch_acc = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0

    print("**********************")
    print("** Training Dataset **")
    print("**********************")
    pbar = tqdm(total=len(train_loader), ncols=0, desc="Train", unit=" step")
    
    A_G_TP = 1e-6
    A_G_FN = 1e-6
    A_G_FP = 1e-6
    A_G_TN = 1e-6
    D_G_TP = 1e-6
    D_G_FN = 1e-6
    D_G_FP = 1e-6
    D_G_TN = 1e-6

    J_G_TP = 1e-6
    J_G_FN = 1e-6
    J_G_FP = 1e-6
    J_G_TN = 1e-6
    #######################################
    # Important => setting model into training mode
    #######################################   
    model.train()
    for batch_idx, data in enumerate(train_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 800, 4])
        # labels:  torch.Size([40, 1, 800, 3])
        DNAs, labels, chr = data 
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model, criterion)


        #######################################
        # predicting all bp.
        #######################################    
        is_expr = (labels.sum(axis=(1,2)) >= 1)
        # print("is_expr: ", is_expr)

        # Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
        A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
        D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
        D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()

        J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, JUNC_THRESHOLD, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        
        A_accuracy, A_auprc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auprc = print_top_1_statistics(Donor_YL, Donor_YP)
        A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, JUNC_THRESHOLD, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
        D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, JUNC_THRESHOLD, D_G_TP, D_G_FN, D_G_FP, D_G_TN)

        batch_loss = loss.item()
        epoch_loss += loss.item()
        epoch_donor_acc += D_accuracy
        epoch_acceptor_acc += A_accuracy

        pbar.update(1)
        pbar.set_postfix(
            batch_id=batch_idx,
            idx_train=len(train_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            # accuracy=f"{batch_acc:.6f}",
            # A_accuracy=f"{A_accuracy:.6f}",
            # D_accuracy=f"{D_accuracy:.6f}",
            A_auprc = f"{A_auprc:.6f}",
            D_auprc = f"{D_auprc:.6f}",
            # A_TP=A_TP,
            # A_FN=A_FN, 
            # A_FP=A_FP, 
            # A_TN=A_TN,
            # D_TP=D_TP,
            # D_FN=D_FN, 
            # D_FP=D_FP, 
            # D_TN=D_TN,
            A_Precision=f"{A_TP/(A_TP+A_FP+1e-6):.6f}",
            A_Recall=f"{A_TP/(A_TP+A_FN+1e-6):.6f}",
            D_Precision=f"{D_TP/(D_TP+D_FP+1e-6):.6f}",
            D_Recall=f"{D_TP/(D_TP+D_FN+1e-6):.6f}",
            J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
            J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
        )
        loss.backward()
        optimizer.step()
        scheduler.step()
        optimizer.zero_grad()

        fw_train_log_loss.write(str(batch_loss)+ "\n")
        fw_train_log_lr.write(str(get_lr(optimizer))+ "\n")
        fw_train_log_A_auprc.write(str(A_auprc)+ "\n")
        fw_train_log_A_threshold_precision.write(f"{A_TP/(A_TP+A_FP+1e-6):.6f}\n")
        fw_train_log_A_threshold_recall.write(f"{A_TP/(A_TP+A_FN+1e-6):.6f}\n")
        fw_train_log_D_auprc.write(str(D_auprc)+ "\n")
        fw_train_log_D_threshold_precision.write(f"{D_TP/(D_TP+D_FP+1e-6):.6f}\n")
        fw_train_log_D_threshold_recall.write(f"{D_TP/(D_TP+D_FN+1e-6):.6f}\n")
        fw_train_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
        fw_train_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
    pbar.close()

    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Donor top-k Acc: {epoch_donor_acc/len(train_loader):.3f} | Acceptor top-k Acc: {epoch_acceptor_acc/len(train_loader):.3f}')
    print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    print(f'Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
    print(f'Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
    print ("Learning rate: %.5f" % (get_lr(optimizer)))
    print("\n\n")


def val_one_epoch(epoch_idx, val_loader):
    epoch_loss = 0
    epoch_acc = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0

    print("**********************")
    print("** Validation Dataset **")
    print("**********************")
    pbar = tqdm(total=len(val_loader), ncols=0, desc="Validation", unit=" step")
    
    A_G_TP = 1e-6
    A_G_FN = 1e-6
    A_G_FP = 1e-6
    A_G_TN = 1e-6
    D_G_TP = 1e-6
    D_G_FN = 1e-6
    D_G_FP = 1e-6
    D_G_TN = 1e-6

    J_G_TP = 1e-6
    J_G_FN = 1e-6
    J_G_FP = 1e-6
    J_G_TN = 1e-6
    #######################################
    # Important => setting model into evaluation mode
    #######################################   
    model.eval()
    for batch_idx, data in enumerate(val_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 800, 4])
        # labels:  torch.Size([40, 1, 800, 3])
        DNAs, labels, chr = data 
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model, criterion)


        #######################################
        # predicting all bp.
        #######################################    
        is_expr = (labels.sum(axis=(1,2)) >= 1)
        # print("is_expr: ", is_expr)

        Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
        A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
        D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
        D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()

        J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, JUNC_THRESHOLD, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        
        A_accuracy, A_auprc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auprc = print_top_1_statistics(Donor_YL, Donor_YP)
        A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, JUNC_THRESHOLD, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
        D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, JUNC_THRESHOLD, D_G_TP, D_G_FN, D_G_FP, D_G_TN)

        batch_loss = loss.item()
        epoch_loss += loss.item()
        epoch_donor_acc += D_accuracy
        epoch_acceptor_acc += A_accuracy

        pbar.update(1)
        pbar.set_postfix(
            batch_id=batch_idx,
            idx_val=len(val_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            # accuracy=f"{batch_acc:.6f}",
            # A_accuracy=f"{A_accuracy:.6f}",
            # D_accuracy=f"{D_accuracy:.6f}",
            A_auprc = f"{A_auprc:.6f}",
            D_auprc = f"{D_auprc:.6f}",
            # A_TP=A_TP,
            # A_FN=A_FN, 
            # A_FP=A_FP, 
            # A_TN=A_TN,
            # D_TP=D_TP,
            # D_FN=D_FN, 
            # D_FP=D_FP, 
            # D_TN=D_TN,
            A_Precision=f"{A_TP/(A_TP+A_FP+1e-6):.6f}",
            A_Recall=f"{A_TP/(A_TP+A_FN+1e-6):.6f}",
            D_Precision=f"{D_TP/(D_TP+D_FP+1e-6):.6f}",
            D_Recall=f"{D_TP/(D_TP+D_FN+1e-6):.6f}",
            J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
            J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
        )

        fw_val_log_loss.write(str(batch_loss)+ "\n")
        fw_val_log_A_auprc.write(str(A_auprc)+ "\n")
        fw_val_log_A_threshold_precision.write(f"{A_TP/(A_TP+A_FP+1e-6):.6f}\n")
        fw_val_log_A_threshold_recall.write(f"{A_TP/(A_TP+A_FN+1e-6):.6f}\n")
        fw_val_log_D_auprc.write(str(D_auprc)+ "\n")
        fw_val_log_D_threshold_precision.write(f"{D_TP/(D_TP+D_FP+1e-6):.6f}\n")
        fw_val_log_D_threshold_recall.write(f"{D_TP/(D_TP+D_FN+1e-6):.6f}\n")
        fw_val_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
        fw_val_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
    pbar.close()

    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(val_loader):.5f} | Donor top-k Acc: {epoch_donor_acc/len(val_loader):.3f} | Acceptor top-k Acc: {epoch_acceptor_acc/len(val_loader):.3f}')
    print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    print(f'Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
    print(f'Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
    print ("Learning rate: %.5f" % (get_lr(optimizer)))
    print("\n\n")





#######################################
# predicting splice / non-splice
#######################################    
# def val_one_epoch(epoch_idx, val_loader):
#     print("************************")
#     print("** Validation Dataset **")
#     print("************************")
#     epoch_loss = 0
#     epoch_acc = 0
#     pbar = tqdm(total=len(val_loader), ncols=0, desc="Validation", unit=" step")

#     J_G_TP = 1e-6
#     J_G_FN = 1e-6
#     J_G_FP = 1e-6
#     J_G_TN = 1e-6
#     #######################################
#     # Important => setting model into evaluation mode
#     #######################################  
#     model.eval()
#     with torch.no_grad():
#         for batch_idx, data in enumerate(val_loader):
#             # print("batch_idx: ", batch_idx)
#             # DNAs:  torch.Size([40, 800, 4])
#             # labels:  torch.Size([40, 1, 800, 3])
#             DNAs, labels, chr = data 
#             DNAs = DNAs.to(torch.float32).to(device)
#             labels = labels.to(torch.float32).to(device)
#             DNAs = torch.permute(DNAs, (0, 2, 1))
#             loss, accuracy, yp = model_fn(DNAs, labels, model, criterion)
            
#             #######################################
#             # predicting splice / non-splice
#             #######################################    
#             batch_loss = loss.item()
#             batch_acc = accuracy
#             epoch_loss += loss.item()
#             epoch_acc += accuracy

#             labels = labels.to("cpu")
#             yp = yp.to("cpu")
#             J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = junc_statistics(labels, yp, JUNC_THRESHOLD, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        

#             pbar.update(1)
#             pbar.set_postfix(
#                 batch_id=batch_idx,
#                 idx_train=len(train_loader)*BATCH_SIZE,
#                 loss=f"{batch_loss:.6f}",
#                 accuracy=f"{batch_acc:.6f}",
#                 J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
#                 J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
#             )

#             fw_train_log_loss.write(str(batch_loss)+ "\n")
#             fw_train_log_acc.write(str(batch_acc)+ "\n")
#             fw_train_log_lr.write(str(get_lr(optimizer))+ "\n")
#             fw_train_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
#             fw_train_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")

#             fw_val_log_loss.write(str(batch_loss)+ "\n")
#             fw_val_log_acc.write(str(batch_acc)+"\n")
#             fw_val_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
#             fw_val_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
#     pbar.close()
#     print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(val_loader):.5f} | Acc: {epoch_acc/len(val_loader):.3f}')
#     print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')







def test_one_epoch(epoch_idx, test_loader):
    epoch_loss = 0
    epoch_acc = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0

    print("**********************")
    print("** Testing Dataset **")
    print("**********************")
    pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")
    
    A_G_TP = 1e-6
    A_G_FN = 1e-6
    A_G_FP = 1e-6
    A_G_TN = 1e-6
    D_G_TP = 1e-6
    D_G_FN = 1e-6
    D_G_FP = 1e-6
    D_G_TN = 1e-6

    J_G_TP = 1e-6
    J_G_FN = 1e-6
    J_G_FP = 1e-6
    J_G_TN = 1e-6
    #######################################
    # Important => setting model into evaluation mode
    #######################################   
    model.eval()
    for batch_idx, data in enumerate(test_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 800, 4])
        # labels:  torch.Size([40, 1, 800, 3])
        DNAs, labels, chr = data 
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model, criterion)


        #######################################
        # predicting all bp.
        #######################################    
        is_expr = (labels.sum(axis=(1,2)) >= 1)
        # print("is_expr: ", is_expr)

        # Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
        A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
        D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
        D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()

        J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, JUNC_THRESHOLD, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        
        A_accuracy, A_auprc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auprc = print_top_1_statistics(Donor_YL, Donor_YP)
        A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, JUNC_THRESHOLD, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
        D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, JUNC_THRESHOLD, D_G_TP, D_G_FN, D_G_FP, D_G_TN)

        batch_loss = loss.item()
        epoch_loss += loss.item()
        epoch_donor_acc += D_accuracy
        epoch_acceptor_acc += A_accuracy

        pbar.update(1)
        pbar.set_postfix(
            batch_id=batch_idx,
            idx_test=len(test_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            # accuracy=f"{batch_acc:.6f}",
            # A_accuracy=f"{A_accuracy:.6f}",
            # D_accuracy=f"{D_accuracy:.6f}",
            A_auprc = f"{A_auprc:.6f}",
            D_auprc = f"{D_auprc:.6f}",
            # A_TP=A_TP,
            # A_FN=A_FN, 
            # A_FP=A_FP, 
            # A_TN=A_TN,
            # D_TP=D_TP,
            # D_FN=D_FN, 
            # D_FP=D_FP, 
            # D_TN=D_TN,
            A_Precision=f"{A_TP/(A_TP+A_FP+1e-6):.6f}",
            A_Recall=f"{A_TP/(A_TP+A_FN+1e-6):.6f}",
            D_Precision=f"{D_TP/(D_TP+D_FP+1e-6):.6f}",
            D_Recall=f"{D_TP/(D_TP+D_FN+1e-6):.6f}",
            J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
            J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
        )

        fw_test_log_loss.write(str(batch_loss)+ "\n")
        fw_test_log_A_auprc.write(str(A_auprc)+ "\n")
        fw_test_log_A_threshold_precision.write(f"{A_TP/(A_TP+A_FP+1e-6):.6f}\n")
        fw_test_log_A_threshold_recall.write(f"{A_TP/(A_TP+A_FN+1e-6):.6f}\n")
        fw_test_log_D_auprc.write(str(D_auprc)+ "\n")
        fw_test_log_D_threshold_precision.write(f"{D_TP/(D_TP+D_FP+1e-6):.6f}\n")
        fw_test_log_D_threshold_recall.write(f"{D_TP/(D_TP+D_FN+1e-6):.6f}\n")
        fw_test_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
        fw_test_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
    pbar.close()

    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Donor top-k Acc: {epoch_donor_acc/len(test_loader):.3f} | Acceptor top-k Acc: {epoch_acceptor_acc/len(test_loader):.3f}')
    print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    print(f'Donor Precision   : {D_G_TP/(D_G_TP+D_G_FP):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
    print(f'Acceptor Precision: {A_G_TP/(A_G_TP+A_G_FP):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')
    print ("Learning rate: %.5f" % (get_lr(optimizer)))
    print("\n\n")
        

def main():
    #############################
    # Model Training
    #############################
    for epoch_num in range(EPOCH_NUM):
        train_one_epoch(epoch_num, train_loader)
        # val_one_epoch(epoch_num, val_loader)
        test_one_epoch(epoch_num, test_loader)
        torch.save(model, MODEL_OUTPUT_BASE+'splam_'+str(epoch_num)+'.pt')

    fw_train_log_loss.close()
    fw_train_log_acc.close()
    fw_train_log_lr.close()

    fw_val_log_loss.close()
    fw_val_log_acc.close()

    fw_test_log_loss.close()
    fw_test_log_acc.close()
    
if __name__ == "__main__":
    main()