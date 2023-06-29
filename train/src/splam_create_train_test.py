from splam_utils import *
from splam_dataset_Chromsome import *
# from splam_Conformer import *
from SPLAM import *
import numpy as np
import torch
import torch.nn as nn
import sys
from torch.optim import Optimizer, AdamW
from torch.optim.lr_scheduler import LambdaLR
from tqdm import tqdm
import math
import os

def main():
    #############################
    # Global variable definition
    #############################
    EPOCH_NUM = 20
    BATCH_SIZE = 100
    N_WORKERS = 1

    # print("\033[1mContext nucleotides: %d\033[0m" % (SL))
    # print("\033[1mSequence length (output): %d\033[0m" % (SL))
    # same_seeds(0)
    # device = torch.device("cuda" if torch.cuda.is_available() else "mps")
    # print(f"[Info]: Use {device} now!")

    L = 64
    W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11,
                    21, 21, 21, 21, 41, 41, 41, 41, 81, 81, 81, 81, 161, 161, 161, 161])
    AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4,
                    10, 10, 10, 10, 25, 25, 25, 25, 50, 50, 50, 50, 100, 100, 100, 100])
    CL = 2 * np.sum(AR*(W-1))
    print("\033[1mContext nucleotides: %d\033[0m" % (SL))
    print("\033[1mSequence length (output): %d\033[0m" % (SL))
    same_seeds(0)
    device = torch.device("cuda" if torch.cuda.is_available() else "mps")
    print(f"[Info]: Use {device} now!")


    #############################
    # Model Initialization
    #############################
    model = SPLAM(L, W, AR).to(device)
    # criterion = nn.CrossEntropyLoss()
    # optimizer = AdamW(model.parameters(), lr=1e-4)
    # scheduler = get_cosine_schedule_with_warmup(optimizer, 1000, 200000)
    print(f"[Info]: Finish creating model!",flush = True)
    print("model: ", model)


    #############################
    # Training Data initialization
    #############################
    save_dataloader(BATCH_SIZE, N_WORKERS)

if __name__ == "__main__":
    main()