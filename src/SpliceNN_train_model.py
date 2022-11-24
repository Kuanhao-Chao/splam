from SpliceNN_utils import *
from SpliceNN_dataset import *
from SpliceNN import *
import numpy as np
import torch
import torch.nn as nn
import sys
from torch.optim import Optimizer, AdamW
from torch.optim.lr_scheduler import LambdaLR
from tqdm import tqdm
import math

#############################
# Global variable definition
#############################
EPOCH_NUM = 10
BATCH_SIZE = 40
N_WORKERS = 1
L = 32

if int(sys.argv[1]) == 80:
    W = np.asarray([11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1])
elif int(sys.argv[1]) == 400:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
elif int(sys.argv[1]) == 2000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                    10, 10, 10, 10])
elif int(sys.argv[1]) == 10000:
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21, 41, 41, 41, 41])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                    10, 10, 10, 10, 25, 25, 25, 25])

CL = 2 * np.sum(AR*(W-1))
assert CL <= CL_MAX and CL == int(sys.argv[1])
print("\033[1mContext nucleotides: %d\033[0m" % (CL))
print("\033[1mSequence length (output): %d\033[0m" % (SL))
same_seeds(0)
device = torch.device("cuda" if torch.cuda.is_available() else "mps")
print(f"[Info]: Use {device} now!")

#############################
# Model Initialization
#############################
model = SpliceNN(L, W, AR).to(device)
criterion = nn.CrossEntropyLoss()
optimizer = AdamW(model.parameters(), lr=1e-4)
# scheduler = get_cosine_schedule_with_warmup(optimizer, 1000, 200000)
print(f"[Info]: Finish creating model!",flush = True)
print("model: ", model)


#############################
# Training Data initialization
#############################
train_loader, test_loader = get_dataloader(BATCH_SIZE, N_WORKERS)
# train_iterator = iter(train_loader)
# valid_iterator = iter(valid_loader)
# print(f"[Info]: Finish loading data!",flush = True)
print("train_iterator: ", len(train_loader))
print("valid_iterator: ", len(test_loader))


def train_one_epoch(epoch_idx, train_loader):
    print("**********************")
    print("** Training Dataset **")
    print("**********************")
    pbar = tqdm(total=len(train_loader), ncols=0, desc="Train", unit=" step")
    for batch_idx, data in enumerate(train_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 800, 4])
        # labels:  torch.Size([40, 1, 800, 3])
        DNAs, labels = data 
        # print("\nDNAs: ", DNAs.size())
        # print("labels: ", labels.size())
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)
        # print("DNAs  : ", DNAs)
        # print("labels: ", labels)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model, criterion, device)
        is_expr = (labels.sum(axis=(1,2)) >= 1)

        Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        batch_loss = loss.item()

        pbar.update(1)
        pbar.set_postfix(
            epoch=batch_idx,
            idx_train=len(train_loader)*BATCH_SIZE,
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
    pbar.close()


def test_one_epoch(epoch_idx, test_loader):
    print("*********************")
    print("** Testing Dataset **")
    print("*********************")

    Y_true_1 = [[] for t in range(1)]
    Y_true_2 = [[] for t in range(1)]
    Y_pred_1 = [[] for t in range(1)]
    Y_pred_2 = [[] for t in range(1)]
    pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")
    for batch_idx, data in enumerate(test_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 800, 4])
        # labels:  torch.Size([40, 1, 800, 3])
        DNAs, labels = data 
        # print("\nDNAs: ", DNAs.size())
        # print("labels: ", labels.size())
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)
        # print("DNAs  : ", DNAs)
        # print("labels: ", labels)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model, criterion, device)
        # Yl.append(labels)
        # Yp.append(yp)
        is_expr = (labels.sum(axis=(1,2)) >= 1)

        Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        batch_loss = loss.item()

        pbar.update(1)
        pbar.set_postfix(
            epoch=batch_idx,
            idx_train=len(train_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            A_accuracy=f"{A_accuracy:.6f}",
            D_accuracy=f"{D_accuracy:.6f}",
            A_auc = f"{A_auc:.6f}",
            D_auc = f"{D_auc:.6f}",
        )
        break
    pbar.close()

    # # Yl = torch.stack(Yl, dim=1)
    # # Yp = torch.stack(Yp, dim=1)
    # print("labels.size(): ", labels.size())
    # print("yp.size(): ", yp.size())

    # for t in range(1):
    #     # print("Yl[t].sum(axis=(1,2)): ", len(Yl[t].sum(axis=(1,2))))
    #     is_expr = (labels.sum(axis=(1,2)) >= 1)
    #     # print("is_expr: ", len(is_expr))
    #     Y_true_1[t].extend(labels[is_expr, 1, :].flatten().to('cpu').detach().numpy())
    #     Y_true_2[t].extend(labels[is_expr, 2, :].flatten().to('cpu').detach().numpy())
    #     Y_pred_1[t].extend(yp[is_expr, 1, :].flatten().to('cpu').detach().numpy())
    #     Y_pred_2[t].extend(yp[is_expr, 2, :].flatten().to('cpu').detach().numpy())

    # print("\n\n##########################################################")
    # # print(">> valid_loader length: ", len(valid_loader))
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
    # print("\n\n")


def main():
    #############################
    # Model Training
    #############################
    for epoch_num in range(EPOCH_NUM):
        train_one_epoch(epoch_num, train_loader)
        test_one_epoch(epoch_num, test_loader)
        torch.save(model, './MODEL/SpliceNN_'+str(epoch_num)+'.pt')
        # test_one_epoch(epoch_num)

if __name__ == "__main__":
    main()