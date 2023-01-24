import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
# from TEST_dataset import *
from SpliceNN_dataset_Chromsome import *
from SpliceNN import *
from SpliceNN_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
from tqdm import tqdm
import warnings
from sklearn.metrics import precision_recall_curve, roc_curve


warnings.filterwarnings("ignore")

#############################
# Global variable definition
#############################
EPOCH_NUM = 20
BATCH_SIZE = 100
N_WORKERS = 1
device = torch.device("cuda" if torch.cuda.is_available() else "mps")
model = torch.load("./MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L600_v20/SpliceNN_19.pt")

#############################
# Model Initialization
#############################
# criterion = nn.CrossEntropyLoss()
print(f"[Info]: Finish loading model!",flush = True)
print("model: ", model)

#############################
# Training Data initialization
#############################
# TARGET = 'positive'
# TARGET = 'negative_canonical'
TARGET = 'test'

# test_loader = get_dataloader(BATCH_SIZE, 'negative_canonical', N_WORKERS)
test_loader = get_dataloader(BATCH_SIZE, N_WORKERS)
# test_loader = get_dataloader(BATCH_SIZE, 'negative_noncanonical', N_WORKERS)

# train_iterator = iter(train_loader)
# valid_iterator = iter(valid_loader)
# print(f"[Info]: Finish loading data!",flush = True)
print("valid_iterator: ", len(test_loader))
MODEL_OUTPUT_BASE = "./TEST/OUTPUT/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L600_v20/"
TARGET_OUTPUT_BASE = MODEL_OUTPUT_BASE + TARGET + "/"
LOG_OUTPUT_TEST_BASE = TARGET_OUTPUT_BASE + "LOG/"
os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)

############################
# Log for testing
############################
test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
test_log_A_topk_accuracy = LOG_OUTPUT_TEST_BASE + "test_A_topk_accuracy.txt"
test_log_A_auc = LOG_OUTPUT_TEST_BASE + "test_A_auc.txt"
test_log_A_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_A_threshold_precision.txt"
test_log_A_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_A_threshold_recall.txt"
test_log_D_topk_accuracy = LOG_OUTPUT_TEST_BASE + "test_D_topk_accuracy.txt"
test_log_D_auc = LOG_OUTPUT_TEST_BASE + "test_D_auc.txt"
test_log_D_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_D_threshold_precision.txt"
test_log_D_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_D_threshold_recall.txt"

fw_test_log_loss = open(test_log_loss, 'w')
fw_test_log_A_topk_accuracy = open(test_log_A_topk_accuracy, 'w')
fw_test_log_A_auc = open(test_log_A_auc, 'w')
fw_test_log_A_threshold_precision = open(test_log_A_threshold_precision, 'w')
fw_test_log_A_threshold_recall = open(test_log_A_threshold_recall, 'w')
fw_test_log_D_topk_accuracy = open(test_log_D_topk_accuracy, 'w')
fw_test_log_D_auc = open(test_log_D_auc, 'w')
fw_test_log_D_threshold_precision = open(test_log_D_threshold_precision, 'w')
fw_test_log_D_threshold_recall = open(test_log_D_threshold_recall, 'w')

def plot_pr_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = precision_recall_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('Recall')
    plt.ylabel('Precision')

def plot_roc_curve(true_y, y_prob, label):
    """
    plots the roc curve based of the probabilities
    """
    fpr, tpr, thresholds = roc_curve(true_y, y_prob)
    plt.plot(fpr, tpr, label=label)
    plt.legend()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')


def test_one_epoch(epoch_idx, test_loader):
    print("*********************")
    print("** Testing Dataset **")
    print("*********************")
    epoch_loss = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0
    pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")

    A_G_TP = 0
    A_G_FN = 0
    A_G_FP = 0
    A_G_TN = 0
    D_G_TP = 0
    D_G_FN = 0
    D_G_FP = 0
    D_G_TN = 0

    J_G_TP = 1e-6
    J_G_FN = 1e-6
    J_G_FP = 1e-6
    J_G_TN = 1e-6
    Acceptor_Sum = np.zeros(600)
    Donor_Sum = np.zeros(600)
    # print("Acceptor_Sum: ", Acceptor_Sum.shape)
    # print("Acceptor_Sum: ", Acceptor_Sum)
    # print("Donor_Sum: ", Donor_Sum.shape)
    # print("Donor_Sum: ", Donor_Sum)
    # Acceptor_Sum = Acceptor_Sum.to(torch.float32)
    # Donor_Sum = Donor_Sum.to(torch.float32)

    # print("Acceptor_Sum: ", Acceptor_Sum.size())
    # print("Donor_Sum   : ", Donor_Sum.size())

    All_Acceptor_YL = []
    All_Acceptor_YP = []
    All_Donor_YL = []
    All_Donor_YP = []

    for batch_idx, data in enumerate(test_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 800, 4])
        # labels:  torch.Size([40, 1, 800, 3])
        # print("len: ", len(data))
        DNAs, labels, seq_names = data 
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model)
        
        is_expr = (labels.sum(axis=(1,2)) >= 1)

        Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].flatten().to('cpu').detach().numpy()

        All_Acceptor_YL.extend(Acceptor_YL)
        All_Acceptor_YP.extend(Acceptor_YP)
        All_Donor_YL.extend(Donor_YL)
        All_Donor_YP.extend(Donor_YP)

        # print("All_Acceptor_YL: ", All_Acceptor_YL)
        # print("All_Acceptor_YL: ", len(All_Acceptor_YL))
        # print("All_Acceptor_YP: ", All_Acceptor_YP)
        # print("All_Acceptor_YP: ", len(All_Acceptor_YP))
        # All_Acceptor_YP.extend(Acceptor_YP)
        # All_Donor_YL.extend(Donor_YL)
        # All_Donor_YP.extend(Donor_YP)
        
        Acceptor_Sum += yp[is_expr, 1, :].sum(axis=0).to('cpu').detach().numpy()
        Donor_Sum += yp[is_expr, 2, :].sum(axis=0).to('cpu').detach().numpy()

        A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
        A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
        D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
        D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()
        # print("Acceptor_Sum: ", Acceptor_Sum.shape)
        # print("Acceptor_Sum: ", Acceptor_Sum)
        # print("Donor_Sum: ", Donor_Sum.shape)
        # print("Donor_Sum: ", Donor_Sum)

        J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = print_junc_statistics(D_YL, A_YL, D_YP, A_YP, 0.5, J_G_TP, J_G_FN, J_G_FP, J_G_TN)
        A_accuracy, A_auc = print_top_1_statistics(Acceptor_YL, Acceptor_YP)
        D_accuracy, D_auc = print_top_1_statistics(Donor_YL, Donor_YP)


        A_G_TP, A_G_FN, A_G_FP, A_G_TN, A_TP, A_FN, A_FP, A_TN = print_threshold_statistics(Acceptor_YL, Acceptor_YP, 0.5, A_G_TP, A_G_FN, A_G_FP, A_G_TN)
        D_G_TP, D_G_FN, D_G_FP, D_G_TN, D_TP, D_FN, D_FP, D_TN = print_threshold_statistics(Donor_YL, Donor_YP, 0.5, D_G_TP, D_G_FN, D_G_FP, D_G_TN)


        batch_loss = loss.item()
        epoch_loss += loss.item()
        epoch_donor_acc += D_accuracy
        epoch_acceptor_acc += A_accuracy

        pbar.update(1)
        pbar.set_postfix(
            epoch=batch_idx,
            idx_train=len(test_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
            # accuracy=f"{batch_acc:.6f}",
            A_accuracy=f"{A_accuracy:.6f}",
            D_accuracy=f"{D_accuracy:.6f}",
            A_auc = f"{A_auc:.6f}",
            D_auc = f"{D_auc:.6f}",
            A_Precision=f"{A_TP/(A_TP+A_FP+1.e-10):.6f}",
            A_Recall=f"{A_TP/(A_TP+A_FN+1.e-10):.6f}",
            D_Precision=f"{D_TP/(D_TP+D_FP+1.e-10):.6f}",
            D_Recall=f"{D_TP/(D_TP+D_FN+1.e-10):.6f}"
        )
        fw_test_log_loss.write(str(batch_loss)+ "\n")
        fw_test_log_A_topk_accuracy.write(str(A_accuracy)+ "\n")
        fw_test_log_A_auc.write(str(A_auc)+ "\n")
        fw_test_log_A_threshold_precision.write(f"{A_TP/(A_TP+A_FP+1.e-10):.6f}\n")
        fw_test_log_A_threshold_recall.write(f"{A_TP/(A_TP+A_FN+1.e-10):.6f}\n")
        fw_test_log_D_topk_accuracy.write(str(D_accuracy)+ "\n")
        fw_test_log_D_auc.write(str(D_auc)+ "\n")
        fw_test_log_D_threshold_precision.write(f"{D_TP/(D_TP+D_FP+1.e-10):.6f}\n")
        fw_test_log_D_threshold_recall.write(f"{D_TP/(D_TP+D_FN+1.e-10):.6f}\n")
    pbar.close()



    plot_roc_curve(All_Acceptor_YL, All_Acceptor_YP, "Acceptor")
    plot_roc_curve(All_Donor_YL, All_Donor_YP, "Donor")
    plt.savefig("output_600bp_roc.png", dpi=300)
    plt.close()

    plot_pr_curve(All_Acceptor_YL, All_Acceptor_YP, "Acceptor")
    plot_pr_curve(All_Donor_YL, All_Donor_YP, "Donor")
    plt.savefig("output_600bp_pr.png", dpi=300)
    plt.close()

    acceptor_scores = Acceptor_Sum / (len(test_loader)*BATCH_SIZE)
    donor_scores = Donor_Sum / (len(test_loader)*BATCH_SIZE)
    # print("acceptor_scores: ", acceptor_scores)
    # print("donor_scores: ", donor_scores)
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Donor Acc: {epoch_donor_acc/len(test_loader):.3f} | Acceptor Acc: {epoch_acceptor_acc/len(test_loader):.3f}')
    print(f'Expected #prediction: {len(test_loader)*BATCH_SIZE+0:03}')
    print(f'Junction Precision  : {J_G_TP/(J_G_TP+J_G_FP):.5f} | Junction Recall: {J_G_TP/(J_G_TP+J_G_FN):.5f} | TP: {J_G_TP} | FN: {J_G_FN} | FP: {J_G_FP} | TN: {J_G_TN}')
    print(f'Donor Precision     : {D_G_TP/(D_G_TP+D_G_FP+1.e-10):.5f} | Donor Recall   : {D_G_TP/(D_G_TP+D_G_FN+1.e-10):.5f} | TP: {D_G_TP} | FN: {D_G_FN} | FP: {D_G_FP} | TN: {D_G_TN}')
    print(f'Acceptor Precision  : {A_G_TP/(A_G_TP+A_G_FP+1.e-10):.5f} | Acceptor Recall: {A_G_TP/(A_G_TP+A_G_FN+1.e-10):.5f} | TP: {A_G_TP} | FN: {A_G_FN} | FP: {A_G_FP} | TN: {A_G_TN}')

    print("")
    print("\n\n")

    y_pos = np.arange(len(Acceptor_Sum))
    fig, ax = plt.subplots(figsize=(12, 6))
    ax2 = ax.twinx()  

    ax.bar(y_pos, acceptor_scores, align='center', alpha=0.5, width=5, color="blue")
    ax2.bar(y_pos, donor_scores, align='center', alpha=0.5, width=5, color="red")

    plt.ylabel('SpliceNN prediction score')
    plt.title('SpliceNN')

    plt.savefig(TARGET_OUTPUT_BASE+"spliceNN_"+TARGET+".png", dpi=300)



def main():
    #############################
    # Model Training
    #############################
    # for epoch_num in range(EPOCH_NUM):
    test_one_epoch(0, test_loader)


if __name__ == "__main__":
    main()