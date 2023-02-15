#! /usr/bin/python3
"""The Dataset and DataLoader class for SPLAM!

    File name: dataset.py
    Author: Kuan-Hao Chao
    Email: kh.chao@cs.jhu.edu
    Date created: 12/20/2022
    Date last modified: 01/14/2023
    Python Version: 3.8
"""

import os, sys
f = sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import torch
from torch.utils.data import Dataset, DataLoader
import pickle
import random
from splam_utils import *

SEQ_LEN = "800"

def split_seq_name(seq):
    return seq[1:]

class myDataset(Dataset):
    """myDataset for SPLAM!
    """
    def __init__(self, type, of, shuffle, segment_len=800):
        self.segment_len = segment_len
        self.data = []
        self.indices = []
        pidx = 0
        with open(of, "r") as f:
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                if pidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif pidx % 2 == 1:
                    seq = line
                    if seq[0] == ">":
                        seq_name = line
                        continue
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    if X.size()[0] != 800:
                        print("The length of input variable 'X' is not 800 (", seq_name, ")")
                        print(X.size())
                    self.data.append([X, Y, seq_name])
                pidx += 1
                if pidx %10000 == 0:
                    print("pidx: ", pidx)

        index_shuf = list(range(len(self.data)))
        if shuffle:
            random.shuffle(index_shuf)
        list_shuf = [self.data[i] for i in index_shuf]
        self.data = list_shuf 
        self.indices = index_shuf
        print("pidx: ", pidx)

    def __len__(self):
        return len(self.data)
 
    def __getitem__(self, index):
        feature = self.data[index][0]
        label = self.data[index][1]
        seq_name = self.data[index][2]
        feature = torch.flatten(feature, start_dim=1)
        return feature, label, seq_name


def get_dataloader(batch_size, n_workers, output_file, shuffle):
    """Generate dataloader"""
    testset = myDataset("test", output_file, shuffle, int(SEQ_LEN))
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        drop_last=True,
        pin_memory=True,
    )
    if batch_size == 1:
        print("shuffle: ", shuffle)
        torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.nobatch.pt")
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.nobatch.indices.pkl", "wb") as f:
            pickle.dump(testset.indices, f)
    elif shuffle:
        print("shuffle: ", shuffle)
        torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.shuffle.pt")
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.shuffle.indices.pkl", "wb") as f:
            pickle.dump(testset.indices, f)
    else:
        print("shuffle: ", shuffle)
        torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.noshuffle.pt")
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.noshuffle.indices.pkl", "wb") as f:
            pickle.dump(testset.indices, f)
    return test_loader


def get_dataloader(batch_size, n_workers, output_file, shuffle, repeat_idx):
    print("output_file: ", output_file)
    testset = myDataset("test", output_file, shuffle, int(SEQ_LEN))
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        drop_last=True,
        pin_memory=True,
    )
    if batch_size == 1:
        print("shuffle: ", shuffle)
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.nobatch.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
            pickle.dump(testset.indices, f)
            
    elif shuffle:
        print("shuffle: ", shuffle)
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.shuffle.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
            pickle.dump(testset.indices, f)

    else:
        print("shuffle: ", shuffle)
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.noshuffle.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
            pickle.dump(testset.indices, f)
    return test_loader