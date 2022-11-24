import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
import json
import random

from SpliceNN_utils import *

class myDataset(Dataset):
    def __init__(self, segment_len=800):
        self.segment_len = segment_len
        self.data = []

        CONSTANT_SIZE = 100000
        #################################
        ## Processing 'POSITIVE' samples
        #################################
        pidx = 0
        with open("./INPUTS/input_pos.shuffle.fa", "r") as f:
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if pidx % 2 == 0:
                    seq_name = line
                elif pidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])
                    # print(X.size())
                    # print(Y.size())
                    self.data.append([X, Y])
                pidx += 1
                if pidx %10000 == 0:
                    print("pidx: ", pidx)
                if pidx > CONSTANT_SIZE:
                    break

        #####################x############
        ## Processing 'NEGATIVE' samples
        #################################
        nidx = 0
        with open("./INPUTS/input_neg.shuffle.fa", "r") as f:
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if nidx % 2 == 0:
                    seq_name = line
                elif nidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '-')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])
                    # print(X.size())
                    # print(Y.size())
                    self.data.append([X, Y])
                nidx += 1
                if nidx %10000 == 0:
                    print("nidx: ", nidx)
                if nidx > CONSTANT_SIZE:
                    break

        #####################x############
        ## Processing 'Non-canonical NEGATIVE' samples
        #################################
        nnidx = 0
        with open("./INPUTS/input_noncan_neg.shuffle.fa", "r") as f:
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if nnidx % 2 == 0:
                    seq_name = line
                elif nnidx % 2 == 1:
                    seq = line
                    # print(seq)
                    X, Y = create_datapoints(seq, '-')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])
                    # print(X.size())
                    # print(Y.size())
                    self.data.append([X, Y])
                nnidx += 1
                if nnidx %10000 == 0:
                    print("nnidx: ", nnidx)
                if nnidx > CONSTANT_SIZE:
                    break
        random.shuffle(self.data)

    def __len__(self):
        return len(self.data)
 
    def __getitem__(self, index):
        # Load preprocessed mel-spectrogram.
        # print("self.data: ", self.data[index])
        feature = self.data[index][0]
        label = self.data[index][1]
        feature = torch.flatten(feature, start_dim=1)
        return feature, label


def get_dataloader(batch_size, n_workers):
    """Generate dataloader"""
    dataset = myDataset(800)
    # Split dataset into training dataset and validation dataset
    trainlen = int(0.9 * len(dataset))
    print("trainlen: ", trainlen)
    lengths = [trainlen, len(dataset) - trainlen]
    trainset, validset = random_split(dataset, lengths)
    # print("trainset: ", trainset)

    train_loader = DataLoader(
        trainset,
        batch_size=batch_size,
        shuffle=True,
        drop_last=True,
        # num_workers=n_workers,
        pin_memory=True,
        # collate_fn=collate_batch,
    )
    test_loader = DataLoader(
        validset,
        batch_size = batch_size,
        # batch_size=len(validset),
        # num_workers=n_workers,
        drop_last=True,
        pin_memory=True,
        # collate_fn=collate_batch,
    )

    return train_loader, test_loader