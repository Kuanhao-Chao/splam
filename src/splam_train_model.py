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
EPOCH_NUM = 25
BATCH_SIZE = 100
N_WORKERS = 1

L = 64

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
# k-fold validation setup
#############################
k=10
splits=KFold(n_splits=k,shuffle=True,random_state=42)

#############################
# Training Data initialization
#############################
save_dataloader(BATCH_SIZE, N_WORKERS)
train_loader, val_loader, test_loader = get_dataloader(BATCH_SIZE, N_WORKERS)
train_iterator = iter(train_loader)
valid_iterator = iter(val_loader)
test_iterator = iter(test_loader)
print(f"\033[1m[Info]: Finish loading data!\033[0m",flush = True)
print("train_iterator: ", len(train_loader))
print("valid_iterator: ", len(val_loader))
print("valid_iterator: ", len(test_loader))



MODEL_OUTPUT_BASE = "./MODEL/SPLAM_v2/"
LOG_OUTPUT_BASE = MODEL_OUTPUT_BASE + "LOG/"
LOG_OUTPUT_TRAIN_BASE = MODEL_OUTPUT_BASE + "LOG/TRAIN/"
LOG_OUTPUT_VAL_BASE = MODEL_OUTPUT_BASE + "LOG/VAL/"
LOG_OUTPUT_TEST_BASE = MODEL_OUTPUT_BASE + "LOG/TEST/"

os.makedirs(LOG_OUTPUT_TRAIN_BASE, exist_ok=True)
os.makedirs(LOG_OUTPUT_VAL_BASE, exist_ok=True)
os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)

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

train_log_J_threshold_precision = LOG_OUTPUT_TRAIN_BASE + "train_J_threshold_precision.txt"
train_log_J_threshold_recall = LOG_OUTPUT_TRAIN_BASE + "train_J_threshold_recall.txt"

fw_train_log_loss = open(train_log_loss, 'w')
fw_train_log_acc = open(train_log_acc, 'w')
fw_train_log_lr = open(train_log_lr, 'w')

fw_train_log_J_threshold_precision = open(train_log_J_threshold_precision, 'w')
fw_train_log_J_threshold_recall = open(train_log_J_threshold_recall, 'w')

############################
# Log for validation
############################
val_log_loss = LOG_OUTPUT_VAL_BASE + "val_loss.txt"
val_log_acc = LOG_OUTPUT_VAL_BASE + "val_accuracy.txt"

val_log_J_threshold_precision = LOG_OUTPUT_VAL_BASE + "val_J_threshold_precision.txt"
val_log_J_threshold_recall = LOG_OUTPUT_VAL_BASE + "val_J_threshold_recall.txt"

fw_val_log_loss = open(val_log_loss, 'w')
fw_val_log_acc = open(val_log_acc, 'w')

fw_val_log_J_threshold_precision = open(val_log_J_threshold_precision, 'w')
fw_val_log_J_threshold_recall = open(val_log_J_threshold_recall, 'w')

############################
# Log for testing
############################
test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
test_log_acc = LOG_OUTPUT_TEST_BASE + "test_accuracy.txt"

test_log_J_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_J_threshold_precision.txt"
test_log_J_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_J_threshold_recall.txt"

fw_test_log_loss = open(test_log_loss, 'w')
fw_test_log_acc = open(test_log_acc, 'w')

fw_test_log_J_threshold_precision = open(test_log_J_threshold_precision, 'w')
fw_test_log_J_threshold_recall = open(test_log_J_threshold_recall, 'w')


def get_lr(optimizer):
    for param_group in optimizer.param_groups:
        return param_group['lr']


def train_one_epoch(epoch_idx, train_loader):
    epoch_loss = 0
    epoch_acc = 0

    print("**********************")
    print("** Training Dataset **")
    print("**********************")
    pbar = tqdm(total=len(train_loader), ncols=0, desc="Train", unit=" step")

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
        loss, accuracy, yp = model_fn(DNAs, labels, model, criterion)
        
        #######################################
        # predicting splice / non-splice
        #######################################    
        batch_loss = loss.item()
        batch_acc = accuracy
        epoch_loss += loss.item()
        epoch_acc += accuracy

        labels = labels.to("cpu")
        yp = yp.to("cpu")
        J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = junc_statistics(labels, yp, 0.5, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        


        pbar.update(1)
        pbar.set_postfix(
            batch_id=batch_idx,
            idx_train=len(train_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            accuracy=f"{batch_acc:.6f}",
            J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
            J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
        )
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()

        fw_train_log_loss.write(str(batch_loss)+ "\n")
        fw_train_log_acc.write(str(batch_acc)+ "\n")
        fw_train_log_lr.write(str(get_lr(optimizer))+ "\n")
        fw_train_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
        fw_train_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
    pbar.close()
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Acc: {epoch_acc/len(train_loader):.3f}')
    print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    print ("Learning rate: %.5f" % (get_lr(optimizer)))
    print("\n\n")


def val_one_epoch(epoch_idx, val_loader):
    print("************************")
    print("** Validation Dataset **")
    print("************************")
    epoch_loss = 0
    epoch_acc = 0
    pbar = tqdm(total=len(val_loader), ncols=0, desc="Validation", unit=" step")

    J_G_TP = 1e-6
    J_G_FN = 1e-6
    J_G_FP = 1e-6
    J_G_TN = 1e-6
    #######################################
    # Important => setting model into evaluation mode
    #######################################  
    model.eval()
    with torch.no_grad():
        for batch_idx, data in enumerate(val_loader):
            # print("batch_idx: ", batch_idx)
            # DNAs:  torch.Size([40, 800, 4])
            # labels:  torch.Size([40, 1, 800, 3])
            DNAs, labels, chr = data 
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)
            DNAs = torch.permute(DNAs, (0, 2, 1))
            loss, accuracy, yp = model_fn(DNAs, labels, model, criterion)
            
            #######################################
            # predicting splice / non-splice
            #######################################    
            batch_loss = loss.item()
            batch_acc = accuracy
            epoch_loss += loss.item()
            epoch_acc += accuracy

            labels = labels.to("cpu")
            yp = yp.to("cpu")
            J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = junc_statistics(labels, yp, 0.5, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        

            pbar.update(1)
            pbar.set_postfix(
                batch_id=batch_idx,
                idx_train=len(train_loader)*BATCH_SIZE,
                loss=f"{batch_loss:.6f}",
                accuracy=f"{batch_acc:.6f}",
                J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
                J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
            )

            fw_train_log_loss.write(str(batch_loss)+ "\n")
            fw_train_log_acc.write(str(batch_acc)+ "\n")
            fw_train_log_lr.write(str(get_lr(optimizer))+ "\n")
            fw_train_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
            fw_train_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")

            fw_val_log_loss.write(str(batch_loss)+ "\n")
            fw_val_log_acc.write(str(batch_acc)+"\n")
            fw_val_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
            fw_val_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
    pbar.close()
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(val_loader):.5f} | Acc: {epoch_acc/len(val_loader):.3f}')
    print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')


def test_one_epoch(epoch_idx, test_loader):
    print("*********************")
    print("** Testing Dataset **")
    print("*********************")
    epoch_loss = 0
    epoch_acc = 0
    pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")

    J_G_TP = 1e-6
    J_G_FN = 1e-6
    J_G_FP = 1e-6
    J_G_TN = 1e-6
    #######################################
    # Important => setting model into evaluation mode
    #######################################  
    model.eval()
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # print("batch_idx: ", batch_idx)
            # DNAs:  torch.Size([40, 800, 4])
            # labels:  torch.Size([40, 1, 800, 3])
            DNAs, labels, chr = data 
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)
            DNAs = torch.permute(DNAs, (0, 2, 1))
            loss, accuracy, yp = model_fn(DNAs, labels, model, criterion)
            
            #######################################
            # predicting splice / non-splice
            #######################################    
            batch_loss = loss.item()
            batch_acc = accuracy
            epoch_loss += loss.item()
            epoch_acc += accuracy

            labels = labels.to("cpu")
            yp = yp.to("cpu")
            J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = junc_statistics(labels, yp, 0.5, J_G_TP, J_G_FN, J_G_FP, J_G_TN)        

            pbar.update(1)
            pbar.set_postfix(
                batch_id=batch_idx,
                idx_train=len(train_loader)*BATCH_SIZE,
                loss=f"{batch_loss:.6f}",
                accuracy=f"{batch_acc:.6f}",
                J_Precision=f"{J_TP/(J_TP+J_FP+1e-6):.6f}",
                J_Recall=f"{J_TP/(J_TP+J_FN+1e-6):.6f}"
            )

            fw_train_log_loss.write(str(batch_loss)+ "\n")
            fw_train_log_acc.write(str(batch_acc)+ "\n")
            fw_train_log_lr.write(str(get_lr(optimizer))+ "\n")
            fw_train_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
            fw_train_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")

            fw_test_log_loss.write(str(batch_loss)+ "\n")
            fw_test_log_acc.write(str(batch_acc)+"\n")
            fw_test_log_J_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1e-6):.6f}\n")
            fw_test_log_J_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1e-6):.6f}\n")
    pbar.close()
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Acc: {epoch_acc/len(test_loader):.3f}')
    print(f'Junction Precision: {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    print("\n\n")


def main():
    #############################
    # Model Training
    #############################
    for epoch_num in range(EPOCH_NUM):
        train_one_epoch(epoch_num, train_loader)
        val_one_epoch(epoch_num, val_loader)
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