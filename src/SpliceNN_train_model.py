from SpliceNN_utils import *
from SpliceNN_dataset import *
from SpliceNN_Conformer import *
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
BATCH_SIZE = 20

print("\033[1mContext nucleotides: %d\033[0m" % (SL))
print("\033[1mSequence length (output): %d\033[0m" % (SL))
same_seeds(0)
device = torch.device("cuda" if torch.cuda.is_available() else "mps")
print(f"[Info]: Use {device} now!")

#############################
# Model Initialization
#############################
print("**train_parse_args()[model_config]: ", train_parse_args()["model_config"]["config3"])
model = SpliceNN_Conformer(train_parse_args()["model_config"]["config3"], d_model=train_parse_args()["model_config"]["config3"]["d_model"]).to(device)
# criterion = nn.CrossEntropyLoss()
# criterion = nn.BCEWithLogitsLoss()
optimizer = AdamW(model.parameters(), lr=1e-4)
scheduler = get_cosine_schedule_with_warmup(optimizer, 1000, 200000)
print(f"[Info]: Finish creating model!",flush = True)
print("model: ", model)


#############################
# Training Data initialization
#############################
train_loader, test_loader = get_dataloader(BATCH_SIZE, train_parse_args()["n_workers"])
# train_iterator = iter(train_loader)
# valid_iterator = iter(valid_loader)
# print(f"[Info]: Finish loading data!",flush = True)
print("train_iterator: ", len(train_loader))
print("valid_iterator: ", len(test_loader))

# train_iterator = iter(train_loader)
# print("train_iterator: ", train_iterator)

def train_one_epoch(epoch_idx, train_loader):
    epoch_loss = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0

    print("**********************")
    print("** Training Dataset **")
    print("**********************")
    pbar = tqdm(total=len(train_loader), ncols=0, desc="Train", unit=" step")

    # for data in train_loader:
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

        # DNAs = torch.permute(DNAs, (0, 2, 1))
        # labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model)
        
        is_expr = (labels.sum(axis=(1,2)) >= 1)
        # print("is_expr: ", is_expr)

        Acceptor_YL = labels[is_expr, :, 1].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, :, 1].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, :, 2].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, :, 2].flatten().to('cpu').detach().numpy()

        # print("Acceptor_YL: ", Acceptor_YL)
        # print("Acceptor_YP: ", Acceptor_YP)

        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        batch_loss = loss.item()
        epoch_loss += loss.item()
        epoch_donor_acc += D_accuracy
        epoch_acceptor_acc += A_accuracy

        # print("batch_loss: ", batch_loss)
        # print("batch_acc: ", batch_acc)

        pbar.update(1)
        pbar.set_postfix(
            epoch=batch_idx,
            idx_train=len(train_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            # accuracy=f"{batch_acc:.6f}",
            A_accuracy=f"{A_accuracy:.6f}",
            D_accuracy=f"{D_accuracy:.6f}",
            A_auc = f"{A_auc:.6f}",
            D_auc = f"{D_auc:.6f}",
        )
        loss.backward()
        optimizer.step()
        scheduler.step()
        optimizer.zero_grad()
    pbar.close()
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Acc: {epoch_donor_acc/len(train_loader):.3f} | Acc: {epoch_acceptor_acc/len(train_loader):.3f}')

    print("\n\n")

    # return epoch_loss, epoch_acc

def test_one_epoch(epoch_idx, test_loader):
    print("*********************")
    print("** Testing Dataset **")
    print("*********************")
    epoch_loss = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0
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

        # DNAs = torch.permute(DNAs, (0, 2, 1))
        # labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model)
        
        is_expr = (labels.sum(axis=(1,2)) >= 1)

        Acceptor_YL = labels[is_expr, :, 1].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, :, 1].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, :, 2].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, :, 2].flatten().to('cpu').detach().numpy()

        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)
        batch_loss = loss.item()
        epoch_loss += loss.item()
        epoch_donor_acc += D_accuracy
        epoch_acceptor_acc += A_accuracy

        # print("batch_loss: ", batch_loss)
        # print("batch_acc: ", batch_acc)

        pbar.update(1)
        pbar.set_postfix(
            epoch=batch_idx,
            idx_train=len(train_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            # accuracy=f"{batch_acc:.6f}",
            A_accuracy=f"{A_accuracy:.6f}",
            D_accuracy=f"{D_accuracy:.6f}",
            A_auc = f"{A_auc:.6f}",
            D_auc = f"{D_auc:.6f}",
        )
    pbar.close()
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(train_loader):.5f} | Acc: {epoch_donor_acc/len(train_loader):.3f} | Acc: {epoch_acceptor_acc/len(train_loader):.3f}')
    print("\n\n")

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
        torch.save(model, './MODEL/Conformer/SpliceNN_'+str(epoch_num)+'.pt')
        # test_one_epoch(epoch_num)

if __name__ == "__main__":
    main()