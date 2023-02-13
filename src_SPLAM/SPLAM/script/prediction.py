# python prediction.py SRR1352129 ../../../src/MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/SpliceNN_19_traced.pt
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import torch.nn as nn
import torch
from TEST_dataset import *
# from splam_dataset_Chromsome_v2 import *
from SpliceNN import *
from splam_utils import *
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
model = torch.load(argv[1])

#############################
# Model Initialization
#############################
print(f"[Info]: Finish loading model!",flush = True)
print("model: ", model)

#############################
# Training Data initialization
#############################
TARGET = 'positive'

input_fasta = ""

test_loader = get_dataloader(BATCH_SIZE, N_WORKERS, "../test_chr9/fasta/junction.fa", True, str(0))

print("valid_iterator: ", len(test_loader))

os.makedirs("../test/output/", exist_ok=True)

MODEL_OUTPUT_BASE = "../test_chr9/OUTPUT/"
TARGET_OUTPUT_BASE = MODEL_OUTPUT_BASE + "/"
LOG_OUTPUT_TEST_BASE = TARGET_OUTPUT_BASE + "LOG/"
os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)

############################
# Log for testing
############################


def test_one_epoch(epoch_idx, test_loader):

    criterion = nn.BCELoss()
    ############################
    # Log for testing
    ############################
    junc_scores = TARGET_OUTPUT_BASE + "junc_scores.bed"
    test_log_loss = LOG_OUTPUT_TEST_BASE + "test_loss.txt"

    fw_junc_scores = open(junc_scores, 'w')
    fw_test_log_loss = open(test_log_loss, 'w')

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
    junc_counter = 0
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            DNAs, labels, seqname = data 
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)

            DNAs = torch.permute(DNAs, (0, 2, 1))
            loss, accuracy, yps = model_fn(DNAs, labels, model, criterion)
        
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


            for idx in range(BATCH_SIZE):
                junction_score = yps[idx]

                chr, start, end, strand = seqname[idx].split(";")
                if strand == "+":
                    fw_junc_scores.write(chr+ "\t"+ start + "\t" + end + "\tJUNC_" + str(junc_counter) + "\t0\t"+ strand+ "\t" + str(junction_score) + "\n")
                elif strand == "-":
                    fw_junc_scores.write(chr+ "\t"+ end + "\t" + start + "\tJUNC_" + str(junc_counter) + "\t0\t"+ strand+ "\t" + str(junction_score) + "\n")
                junc_counter += 1
    pbar.close()

    fw_junc_scores.close()
    fw_test_log_loss.close()


    print(f'Expected #prediction: {len(test_loader)*BATCH_SIZE+0:03}')

    print("")
    print("\n\n")

def main():
    #############################
    # Model Training
    #############################
    # for epoch_num in range(EPOCH_NUM):
    test_one_epoch(0, test_loader)


if __name__ == "__main__":
    main()