import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
import json
import random
import os
import math

from SpliceNN_utils import *

# Random_90_10 / Chromosome_90_10
# TARGET = "Chromosome_split_p_n_nn_n1"
SEQ_LEN = "800"
os.makedirs("./INPUTS/SPLAM/", exist_ok=True)

def split_seq_name(seq):
    return seq[1:]

class myDataset(Dataset):
    def __init__(self, type, segment_len=800):
        self.segment_len = segment_len
        self.data = []

        if type == "train":
            CONSTANT_SIZE = 176086
        else:
            CONSTANT_SIZE = 23914

        CONSTANT_SIZE_NEG = math.ceil(CONSTANT_SIZE*2/3)
        #################################
        ## Processing 'POSITIVE' samples
        #################################
        pidx = 0
        with open("./OUTPUT/splam.juncs.seq.fa", "r") as f:
            print("./OUTPUT/splam.juncs.seq.fa")
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
                    if pidx < 6000 :
                        X, Y = create_datapoints(seq, '+')
                    else:
                        X, Y = create_datapoints(seq, '-')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])                        
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
                if pidx > CONSTANT_SIZE:
                    break

    def __len__(self):
        return len(self.data)
 
    def __getitem__(self, index):
        # Load preprocessed mel-spectrogram.
        # print("self.data: ", self.data[index])
        feature = self.data[index][0]
        label = self.data[index][1]
        seq_name = self.data[index][2]
        feature = torch.flatten(feature, start_dim=1)
        return feature, label, seq_name

def get_dataloader(batch_size, n_workers):
    """Generate dataloader"""
    testset = myDataset("test", int(SEQ_LEN))

    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        # batch_size=len(validset),
        # num_workers=n_workers,
        drop_last=True,
        pin_memory=True,
        # collate_fn=collate_batch,
    )
    torch.save(test_loader, "./INPUTS/SPLAM/test.pt")
    return test_loader

# def get_dataloader(batch_size, n_workers):
#     # print("Loading dataset: ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
#     train_loader = torch.load("./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")

#     print("Loading dataset: ", "./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
#     test_loader = torch.load("./INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
#     # return test_loader
#     return train_loader, test_loader