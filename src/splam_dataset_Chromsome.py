import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
import json
import random
import os
import math

from splam_utils import *

TARGET = "DataLoader_p_n_nn"
SEQ_LEN = "800"
os.makedirs("./INPUTS/"+SEQ_LEN+"bp/"+TARGET, exist_ok=True)

def split_seq_name(seq):
    return seq[1:]

class myDataset(Dataset):
    def __init__(self, type, segment_len=800):
        self.segment_len = segment_len
        self.data = []

        #################################
        ## Processing 'POSITIVE' samples
        #################################
        pidx = 0
        with open("./INPUTS/"+SEQ_LEN+"bp/input_pos/"+type+"_pos.shuffle.fa", "r") as f:
            print("Processing ./INPUTS/"+SEQ_LEN+"bp/input_pos/"+type+"_pos.shuffle.fa")
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if pidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif pidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    # Y = torch.Tensor(np.array(Y)[0])
                    # print(X)
                    # print(Y)
                    if X.size()[0] != 800:
                        print("seq_name: ", seq_name)
                        print(X.size())
                        print(Y.size())
                    self.data.append([X, Y, seq_name])
                pidx += 1
                if pidx %10000 == 0:
                    print("pidx: ", pidx)
                    # print(seq_name)
                # if pidx > CONSTANT_SIZE:
                #     break
                # if pidx >= 1000:
                #     break
        print("pidx: ", pidx)

        CONSTANT_SIZE = pidx
        # CONSTANT_SIZE_NEG = math.ceil(CONSTANT_SIZE*2/3)
        CONSTANT_SIZE_NEG = CONSTANT_SIZE
        print("\033[1m[INFO] CONSTANT_SIZE     : ", CONSTANT_SIZE, "\033[0m")
        print("\033[1m[INFO] CONSTANT_SIZE_NEG : ", CONSTANT_SIZE_NEG, "\033[0m")

        # #################################
        # ## Processing 'NEGATIVE_1' samples
        # #################################
        # n1idx = 0
        # with open("./INPUTS/"+SEQ_LEN+"bp/input_neg_1/"+type+"_neg_1.shuffle.fa", "r") as f:
        #     print("Processing ./INPUTS/"+SEQ_LEN+"bp/input_neg_1/"+type+"_neg_1.shuffle.fa")
        #     lines = f.read().splitlines()
        #     seq_name = ""
        #     seq = ""
        #     for line in lines:
        #         # print(line)
        #         if n1idx % 2 == 0:
        #             seq_name = split_seq_name(line)
        #         elif n1idx % 2 == 1:
        #             seq = line
        #             # print(seq)
        #             X, Y = create_datapoints(seq, '-')
        #             X = torch.Tensor(np.array(X))
        #             # Y = torch.Tensor(np.array(Y)[0])
        #             # print(X)
        #             # print(Y)
        #             if X.size()[0] != 800:
        #                 print("seq_name: ", seq_name)
        #                 print(X.size())
        #                 print(Y.size())
        #             self.data.append([X, Y, seq_name])
        #         n1idx += 1
        #         if n1idx %10000 == 0:
        #             print("n1idx: ", n1idx)
        #             print(seq_name)
        #         if n1idx >= CONSTANT_SIZE_NEG:
        #             break

        #################################
        ## Processing 'NEGATIVE' samples
        #################################
        nidx = 0
        with open("./INPUTS/"+SEQ_LEN+"bp/input_neg_can/"+type+"_neg_can.shuffle.fa", "r") as f:
            print("Processing ./INPUTS/"+SEQ_LEN+"bp/input_neg_can/"+type+"_neg_can.shuffle.fa")
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if nidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif nidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '-')
                    X = torch.Tensor(np.array(X))
                    # Y = torch.Tensor(np.array(Y)[0])
                    # print(X)
                    # print(Y)
                    if X.size()[0] != 800:
                        print("seq_name: ", seq_name)
                        print(X.size())
                        print(Y.size())
                    self.data.append([X, Y, seq_name])
                nidx += 1
                if nidx %10000 == 0:
                    print("nidx: ", nidx)
                    # print(seq_name)
                if nidx >= CONSTANT_SIZE_NEG:
                    break
        print("nidx: ", nidx)

        #################################
        ## Processing 'Non-canonical NEGATIVE' samples
        #################################
        nnidx = 0
        with open("./INPUTS/"+SEQ_LEN+"bp/input_neg_noncan/"+type+"_neg_noncan.shuffle.fa", "r") as f:
            print("Processing ./INPUTS/"+SEQ_LEN+"bp/input_neg_noncan/"+type+"_neg_noncan.shuffle.fa")
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if nnidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif nnidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '-')
                    X = torch.Tensor(np.array(X))
                    # Y = torch.Tensor(np.array(Y)[0])
                    # print(X)
                    # print(Y)
                    if X.size()[0] != 800:
                        print("seq_name: ", seq_name)
                        print(X.size())
                        print(Y.size())
                    self.data.append([X, Y, seq_name])
                nnidx += 1
                if nnidx %10000 == 0:
                    print("nnidx: ", nnidx)
                    # print(seq_name)
                if nnidx >= CONSTANT_SIZE_NEG:
                    break
        print("nnidx: ", nnidx)
        random.shuffle(self.data)

    def __len__(self):
        return len(self.data)
 
    def __getitem__(self, index):
        # Load preprocessed mel-spectrogram.
        feature = self.data[index][0]
        label = self.data[index][1]
        seq_name = self.data[index][2]
        feature = torch.flatten(feature, start_dim=1)
        return feature, label, seq_name


def save_dataloader(batch_size, n_workers):
    """Generate dataloader"""
    trainset_origin = myDataset("train", int(SEQ_LEN))
    trainset, valset = torch.utils.data.random_split(trainset_origin, [0.9, 0.1])
    testset = myDataset("test", int(SEQ_LEN))
    train_loader = DataLoader(
        trainset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
        pin_memory=True,
    )
    val_loader = DataLoader(
        valset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=False,
        pin_memory=True,
    )
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        drop_last=False,
        pin_memory=True,
    )

    #######################################
    # predicting splice / non-splice
    #######################################
    torch.save(train_loader, "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
    torch.save(val_loader, "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/val.pt")
    torch.save(test_loader, "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")


def get_dataloader(batch_size, n_workers):
    #######################################
    # predicting splice / non-splice
    #######################################
    print("[INFO] Loading dataset: ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
    train_loader = torch.load("./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
    print("[INFO] Loading dataset: ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/val.pt")
    val_loader = torch.load("./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/val.pt")
    print("[INFO] Loading dataset: ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
    test_loader = torch.load("./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
    return train_loader, val_loader, test_loader
