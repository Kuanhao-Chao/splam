#! /usr/bin/python3
"""All the utility function for SPLAM!

    File name: prediction.py
    Author: Kuan-Hao Chao
    Email: kh.chao@cs.jhu.edu
    Date created: 12/20/2022
    Date last modified: 01/14/2023
    Python Version: 3.8
"""

import os, sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))
import torch
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
from progress.bar import Bar
import warnings
from splam.dataset import *
from SpliceNN import *
from splam.splam_utils import *

warnings.filterwarnings("ignore")
argv = sys.argv[1:]

print("")
print("################################")
print("## Start the predictions now! ##")
print("################################")

#############################
# Global variable definition
#############################
BATCH_SIZE = 100
N_WORKERS = None

device = torch.device("cuda" if torch.cuda.is_available() else "mps")
model = torch.load("../../src/MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/SpliceNN_19.pt")
# model = torch.load(argv[1])

#############################
# Model Initialization
#############################
print(f"[Info]: Finish loading model!",flush = True)
print("SPLAM! model: ", model)

#############################
# Training Data initialization
#############################
test_loader = get_dataloader(BATCH_SIZE, N_WORKERS, "./test_chr9/fasta/junction.fa", True, str(0))

MODEL_OUTPUT_BASE = "../test_chr9/OUTPUT/"
TARGET_OUTPUT_BASE = MODEL_OUTPUT_BASE + "/"
LOG_OUTPUT_TEST_BASE = TARGET_OUTPUT_BASE + "LOG/"
os.makedirs(LOG_OUTPUT_TEST_BASE, exist_ok=True)


def test_model(epoch_idx, test_loader):
    criterion = torch.nn.BCELoss()
    ############################
    # Log for testing
    ############################
    junc_scores = TARGET_OUTPUT_BASE + "junc_scores.bed"
    fw_junc_scores = open(junc_scores, 'w')
    print("*********************")
    print("** Testing Dataset **")
    print("*********************")
    model.eval()
    junc_counter = 0    
    pbar = Bar('SPLAM! prediction', max=len(test_loader))
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            DNAs, labels, seqname = data 
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)
            DNAs = torch.permute(DNAs, (0, 2, 1))
            loss, accuracy, yps = model_fn(DNAs, labels, model, criterion)
            labels = labels.to("cpu").detach().numpy()
            yps = yps.to("cpu").detach().numpy()
            pbar.next()

            for idx in range(BATCH_SIZE):
                junction_score = yps[idx]
                chr, start, end, strand = seqname[idx].split(";")
                if strand == "+":
                    fw_junc_scores.write(chr+ "\t"+ start + "\t" + end + "\tJUNC_" + str(junc_counter) + "\t0\t"+ strand+ "\t" + str(junction_score) + "\n")
                elif strand == "-":
                    fw_junc_scores.write(chr+ "\t"+ end + "\t" + start + "\tJUNC_" + str(junc_counter) + "\t0\t"+ strand+ "\t" + str(junction_score) + "\n")
                junc_counter += 1

    pbar.finish()
    fw_junc_scores.close()
    print(f'Expected #prediction: {len(test_loader)*BATCH_SIZE+0:03}')
    print("\n\n\n")


if __name__ == "__main__":
    test_model(0, test_loader)

# python prediction.py SRR1352129 ../../../src/MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/SpliceNN_19_traced.pt
