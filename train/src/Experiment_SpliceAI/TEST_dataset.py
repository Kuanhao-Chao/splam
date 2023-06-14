import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
import json
import random

from SpliceNN_utils import *

class myDataset(Dataset):
    def __init__(self, targe, fa_file, segment_len=600):
        self.segment_len = segment_len
        self.data = []

        # CONSTANT_SIZE = "all"
        CONSTANT_SIZE = "all"
        #################################
        ## Processing 'POSITIVE' samples
        #################################
        idx = 0
        with open(fa_file, "r") as f:
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if idx % 2 == 0:
                    seq_name = line
                    print("seq_name: ", seq_name)
                elif idx % 2 == 1:
                    seq = line
                    if seq[0] == ">":
                        seq_name = line
                        continue
                    # print(seq)
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])
                    # print(X)
                    # print(Y)
                    # print(X.size())
                    # print(Y.size())
                    self.data.append([X, Y, seq_name])
                idx += 1
                if idx %10000 == 0:
                    print("idx: ", idx)
                if CONSTANT_SIZE != "all":
                    if idx > CONSTANT_SIZE:
                        break
        # random.shuffle(self.data)

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


def get_dataloader(batch_size, target, fa_file, n_workers):
    """Generate dataloader"""
    dataset = myDataset(target, fa_file, 600)

    test_loader = DataLoader(
        dataset,
        batch_size = batch_size,
        # batch_size=len(validset),
        # num_workers=n_workers,
        drop_last=True,
        pin_memory=True,
        # collate_fn=collate_batch,
    )

    return test_loader
    