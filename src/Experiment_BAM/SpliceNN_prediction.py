import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch
from TEST_dataset import *
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
SEQ_LEN="1000"

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
test_loader = get_dataloader(BATCH_SIZE, TARGET, "../../results/"+SEQ_LEN+"bp/"+argv[0]+"/INPUTS/input.fa", N_WORKERS)
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

    Acceptor_Sum = np.zeros(1000)
    Donor_Sum = np.zeros(1000)

    threshold = 0.3
    num_good_juncs = 0
    num_bad_juncs = 0

    junc_counter = 0
    for batch_idx, data in enumerate(test_loader):
        # print("batch_idx: ", batch_idx)
        # DNAs:  torch.Size([40, 1000, 4])
        # labels:  torch.Size([40, 1, 1000, 3])
        DNAs, labels, seq_names = data 
        DNAs = DNAs.to(torch.float32).to(device)
        labels = labels.to(torch.float32).to(device)

        DNAs = torch.permute(DNAs, (0, 2, 1))
        labels = torch.permute(labels, (0, 2, 1))
        loss, yp = model_fn(DNAs, labels, model)
        
        is_expr = (labels.sum(axis=(1,2)) >= 1)

        # Acceptor_YL = labels[is_expr, 1, :].flatten().to('cpu').detach().numpy()
        Acceptor_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
        # Donor_YL = labels[is_expr, 2, :].flatten().to('cpu').detach().numpy()
        Donor_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()

        for idx in range(BATCH_SIZE):
            d_idx = [i for i in range(len(Donor_YP[idx])) if Donor_YP[idx][i] > threshold]
            a_idx = [i for i in range(len(Acceptor_YP[idx])) if Acceptor_YP[idx][i] > threshold]

            donor_score = Donor_YP[idx][250]
            acceptor_score = Acceptor_YP[idx][750]

            chr, start, end, strand = seq_names[idx].split(";")
            if 250 in d_idx and 750 in a_idx:
                # print("d_idx: ", d_idx)
                # print("a_idx: ", a_idx)
                num_good_juncs += 1
            else:
                num_bad_juncs += 1
                if strand == "+":
                    fw_removed_juncs.write(chr[1:]+ "\t"+ start+ "\t"+ end+ "\tJUNC\t0\t"+ strand+ "\n")
                elif strand == "-":
                    fw_removed_juncs.write(chr[1:]+ "\t"+ end + "\t" + start + "\tJUNC\t0\t"+ strand+ "\n")
                # print("seq_names: ", chr[1:], start, end, strand)
            

            if strand == "+":
                fw_junc_scores.write(chr[1:]+ "\t"+ start + "\t" + end + "\tJUNC_" + str(junc_counter) + "\t0\t"+ strand+ "\t" + str(donor_score) + "\t" + str(acceptor_score) + "\n")
            elif strand == "-":
                fw_junc_scores.write(chr[1:]+ "\t"+ end + "\t" + start + "\tJUNC_" + str(junc_counter) + "\t0\t"+ strand+ "\t" + str(donor_score) + "\t" + str(acceptor_score) + "\n")
            junc_counter += 1

        Acceptor_Sum += yp[is_expr, 1, :].sum(axis=0).to('cpu').detach().numpy()
        Donor_Sum += yp[is_expr, 2, :].sum(axis=0).to('cpu').detach().numpy()

        # print("Acceptor_Sum: ", Acceptor_Sum.shape)
        # print("Acceptor_Sum: ", Acceptor_Sum)
        # print("Donor_Sum: ", Donor_Sum.shape)
        # print("Donor_Sum: ", Donor_Sum)

        batch_loss = loss.item()
        epoch_loss += loss.item()

        pbar.update(1)
        pbar.set_postfix(
            epoch=batch_idx,
            idx_train=len(test_loader)*BATCH_SIZE,
            loss=f"{batch_loss:.6f}",
        )
        fw_test_log_loss.write(str(batch_loss)+ "\n")
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