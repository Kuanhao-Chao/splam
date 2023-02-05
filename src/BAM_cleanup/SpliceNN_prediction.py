import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch.nn as nn
import torch
# from TEST_dataset import *
from SpliceNN_dataset_Chromsome_v2 import *
from SpliceNN import *
from SpliceNN_utils import *
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
from tqdm import tqdm
import warnings

warnings.filterwarnings("ignore")
argv = sys.argv[1:]
#############################
# Global variable definition
#############################
EPOCH_NUM = 20
BATCH_SIZE = 100
N_WORKERS = 1
SEQ_LEN="800"
QUATER_SEQ_LEN = int(SEQ_LEN)//4

device = torch.device("cuda" if torch.cuda.is_available() else "mps")
model = torch.load("../MODEL/"+argv[1]+"/SpliceNN_19.pt")

#############################
# Model Initialization
#############################
print(f"[Info]: Finish loading model!",flush = True)
print("model: ", model)

#############################
# Training Data initialization
#############################
TARGET = 'positive'

# test_loader = get_dataloader(BATCH_SIZE, 'negative_canonical', N_WORKERS)
# test_loader = get_dataloader(BATCH_SIZE, TARGET, "../../results/"+SEQ_LEN+"bp/"+argv[0]+"/INPUTS/input.fa", N_WORKERS)

test_loader = get_dataloader(BATCH_SIZE, N_WORKERS, "../../results/"+SEQ_LEN+"bp/"+argv[0]+"/INPUTS/input.fa", True, str(0))


# test_loader = get_dataloader(BATCH_SIZE, 'negative_noncanonical', N_WORKERS)

# train_iterator = iter(train_loader)
# valid_iterator = iter(valid_loader)
# print(f"[Info]: Finish loading data!",flush = True)
print("valid_iterator: ", len(test_loader))
MODEL_OUTPUT_BASE = "../../results/"+SEQ_LEN+"bp/"+argv[0]+"/OUTPUT/"+argv[1]+"/"
TARGET_OUTPUT_BASE = MODEL_OUTPUT_BASE + "/"
LOG_OUTPUT_TEST_BASE = TARGET_OUTPUT_BASE + "LOG/"
os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)

############################
# Log for testing
############################
test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
removed_juncs = TARGET_OUTPUT_BASE + "removed_junc.bed"
junc_scores = TARGET_OUTPUT_BASE + "junc_scores.bed"

fw_test_log_loss = open(test_log_loss, 'w')
fw_removed_juncs = open(removed_juncs, 'w')
fw_junc_scores = open(junc_scores, 'w')

def test_one_epoch(epoch_idx, test_loader):
    print("*********************")
    print("** Testing Dataset **")
    print("*********************")
    epoch_loss = 0
    epoch_donor_acc = 0
    epoch_acceptor_acc = 0
    pbar = tqdm(total=len(test_loader), ncols=0, desc="Test", unit=" step")

    Acceptor_Sum = np.zeros(int(SEQ_LEN))
    Donor_Sum = np.zeros(int(SEQ_LEN))

    threshold = 0.3
    num_good_juncs = 0
    num_bad_juncs = 0

    junc_counter = 0
    criterion = nn.BCELoss()

    ############################
    # Log for testing
    ############################
    test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"
    test_log_D_threshold_accuracy = LOG_OUTPUT_TEST_BASE + "test_D_threshold_accuracy.txt"
    test_log_D_threshold_precision = LOG_OUTPUT_TEST_BASE + "test_D_threshold_precision.txt"
    test_log_D_threshold_recall = LOG_OUTPUT_TEST_BASE + "test_D_threshold_recall.txt"

    fw_test_log_loss = open(test_log_loss, 'w')
    fw_test_log_D_threshold_accuracy = open(test_log_D_threshold_accuracy, 'w')
    fw_test_log_D_threshold_precision = open(test_log_D_threshold_precision, 'w')
    fw_test_log_D_threshold_recall = open(test_log_D_threshold_recall, 'w')
    # test_one_epoch(0, test_loader)
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
    All_Junction_YL = []
    All_Junction_YP = []

    model.eval()
    # model.train()
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # print("batch_idx: ", batch_idx)
            # DNAs:  torch.Size([40, 800, 4])
            # labels:  torch.Size([40, 1, 800, 3])
            # print("len: ", len(data))
            DNAs, labels, chr = data 
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)

            DNAs = torch.permute(DNAs, (0, 2, 1))
            # labels = torch.permute(labels, (0, 2, 1))
            loss, accuracy, yps = model_fn(DNAs, labels, model, criterion)
            # loss, accuracy, yps = model_fn(DNAs, labels, model, criterion)
        

            #######################################
            # predicting splice / non-splice
            #######################################    
            batch_loss = loss.item()
            batch_acc = accuracy
            epoch_loss += loss.item()
            epoch_acc += accuracy

            labels = labels.to("cpu").detach().numpy()
            yps = yps.to("cpu").detach().numpy()

            J_G_TP, J_G_FN, J_G_FP, J_G_TN, J_TP, J_FN, J_FP, J_TN = junc_statistics(labels, yps, 0.5, J_G_TP, J_G_FN, J_G_FP, J_G_TN)      

            # labels = [label.to("cpu").detach() for label in labels]
            # yps = [yp.to("cpu").detach() for yp in yps]

            All_Junction_YL.extend(labels)
            All_Junction_YP.extend(yps)


            pbar.update(1)
            pbar.set_postfix(
                epoch=batch_idx,
                idx_train=len(test_loader)*BATCH_SIZE,
                loss=f"{batch_loss:.6f}",
                accuracy=f"{batch_acc:.6f}",
                J_Precision=f"{J_TP/(J_TP+J_FP+1.e-10):.6f}",
                J_Recall=f"{J_TP/(J_TP+J_FN+1.e-10):.6f}"
            )
            fw_test_log_loss.write(str(batch_loss)+ "\n")
            fw_test_log_D_threshold_accuracy.write(str(batch_acc)+ "\n")
            fw_test_log_D_threshold_precision.write(f"{J_TP/(J_TP+J_FP+1.e-10):.6f}\n")
            fw_test_log_D_threshold_recall.write(f"{J_TP/(J_TP+J_FN+1.e-10):.6f}\n")
    pbar.close()














    fw_test_log_loss.close()
    fw_removed_juncs.close()

    acceptor_scores = Acceptor_Sum / (len(test_loader)*BATCH_SIZE)
    donor_scores = Donor_Sum / (len(test_loader)*BATCH_SIZE)
    # print("acceptor_scores: ", acceptor_scores)
    # print("donor_scores: ", donor_scores)
    print(f'Epoch {epoch_idx+0:03}: | Loss: {epoch_loss/len(test_loader):.5f} | Donor Acc: {epoch_donor_acc/len(test_loader):.3f} | Acceptor Acc: {epoch_acceptor_acc/len(test_loader):.3f}')
    print(f'Expected #prediction: {len(test_loader)*BATCH_SIZE+0:03}')
    print("Number of good junctions : ", num_good_juncs)
    print("Number of bad junctions  : ", num_bad_juncs)

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