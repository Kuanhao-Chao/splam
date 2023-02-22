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
import argparse
import numpy as np
import warnings
import matplotlib.pyplot as plt; plt.rcdefaults()
from progress.bar import Bar
from splam.dataset import *
from SPLAM import *
from splam.splam_utils import *

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='SPLAM! splice junction prediction.')
parser.add_argument('-f', metavar='<junction.fa>', required=True, help='the junction FASTA file in SPLAM! format')
parser.add_argument('-o', metavar='<score.bed>', required=True, help='the output SPLAM! scores for junctions')
parser.add_argument('-m', metavar='<model.pt>', required=True, help='the path to the SPLAM! model')

args = parser.parse_args()
JUNC_FA = args.f
OUT_SCORE = args.o
MODEL_PATH = args.m

def test_model():
    #############################
    # Global variable definition
    #############################
    BATCH_SIZE = 100
    N_WORKERS = None

    device = torch.device("cuda" if torch.cuda.is_available() else "mps")
    model = torch.load(MODEL_PATH)

    #############################
    # Model Initialization
    #############################
    print(f"[Info]: Finish loading model!",flush = True)
    print(f"[Info]: SPLAM! prediction",flush = True)
    # print("SPLAM! model: ", model)

    #############################
    # Training Data initialization
    #############################
    # "./test_chr9/fasta/junction.fa"
    test_loader = get_dataloader(BATCH_SIZE, N_WORKERS, JUNC_FA, True, str(0))

    # MODEL_OUTPUT_BASE = "../test_chr9/OUTPUT/"
    # TARGET_OUTPUT_BASE = MODEL_OUTPUT_BASE + "/"

    criterion = torch.nn.BCELoss()

    fw_junc_scores = open(OUT_SCORE, 'w')
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
            for idx in range(len(yps)):
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
    print("\n\n")


if __name__ == "__main__":
    test_model()

# python prediction.py SRR1352129 ../../../src/MODEL/SpliceAI_6_RB_p_n_nn_n1_TB_all_samples_thr_100_splitByChrom_L64_C16_L800_v31/SpliceNN_19_traced.pt
