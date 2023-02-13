import torch
from torch.utils.data import Dataset, DataLoader, random_split
from pathlib import Path
import json
import random
import os
import math
import pickle

from splam_utils import *

# Random_90_10 / Chromosome_90_10
# TARGET = "Chromosome_split_p_n_nn_n1"
SEQ_LEN = "800"
os.makedirs("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/", exist_ok=True)

def split_seq_name(seq):
    return seq[1:]



class myDataset(Dataset):
    def __init__(self, type, of, shuffle, segment_len=800):
        self.segment_len = segment_len
        self.data = []
        self.indices = []
        if type == "train":
            CONSTANT_SIZE = 176086
        else:
            CONSTANT_SIZE = 23914

        # CONSTANT_SIZE = 500
        CONSTANT_SIZE_NEG = math.ceil(CONSTANT_SIZE*2/3)
        #################################
        ## Processing 'POSITIVE' samples
        #################################
        pidx = 0


        with open(of, "r") as f:
            print("of: ", of)
            lines = f.read().splitlines()
            seq_name = ""
            seq = ""
            for line in lines:
                # print(line)
                if pidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif pidx % 2 == 1:
                    seq = line
                    if seq[0] == ">":
                        seq_name = line
                        continue
                    
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    # print(X)
                    # print(Y)
                    if X.size()[0] != 800:
                        print("seq_name: ", seq_name)
                        print(X.size())
                        # print(Y.size())
                    self.data.append([X, Y, seq_name])
                pidx += 1
                if pidx %10000 == 0:
                    print("pidx: ", pidx)

        index_shuf = list(range(len(self.data)))
        if shuffle:
            random.shuffle(index_shuf)
            # Shuffle just in a certain range.

        list_shuf = [self.data[i] for i in index_shuf]
        self.data = list_shuf 
        self.indices = index_shuf
        print("pidx: ", pidx)

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



def get_dataloader(batch_size, n_workers, output_file, shuffle):
    """Generate dataloader"""
    # testset = myDataset("test", output_file, int(SEQ_LEN))
    testset = myDataset("test", output_file, shuffle, int(SEQ_LEN))

    print("testset.indices: ", testset.indices)

    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        # batch_size=len(validset),
        # num_workers=n_workers,
        drop_last=True,
        pin_memory=True,
        # collate_fn=collate_batch,
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

# def get_dataloader(batch_size, n_workers):
#     # print("Loading dataset: ", "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")
#     train_loader = torch.load("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/train.pt")

#     print("Loading dataset: ", "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
#     test_loader = torch.load("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/"+SEQ_LEN+"bp/"+TARGET+"/test.pt")
#     # return test_loader
#     return train_loader, test_loader

def get_dataloader(batch_size, n_workers, output_file, shuffle, repeat_idx):
    """Generate dataloader"""
    # testset = myDataset("test", output_file, int(SEQ_LEN))
    print("output_file: ", output_file)
    testset = myDataset("test", output_file, shuffle, int(SEQ_LEN))

    # print("testset.indices: ", testset.indices)

    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        # batch_size=len(validset),
        # num_workers=n_workers,
        drop_last=True,
        pin_memory=True,
        # collate_fn=collate_batch,
    )
    if batch_size == 1:
        print("shuffle: ", shuffle)
        # torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.nobatch."+str(repeat_idx)+"."+str(batch_size)+".pt")
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.nobatch.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
            pickle.dump(testset.indices, f)
            
    elif shuffle:
        print("shuffle: ", shuffle)
        # torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.shuffle."+str(repeat_idx)+"."+str(batch_size)+".pt")
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.shuffle.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
            pickle.dump(testset.indices, f)

    else:
        print("shuffle: ", shuffle)
        # torch.save(test_loader, "/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUTS/SPLAM_v2/test.noshuffle."+str(repeat_idx)+"."+str(batch_size)+".pt")
        with open("/Users/chaokuan-hao/Documents/Projects/PR_SPLAM/src_tools_evaluation/INPUT/splam.noshuffle.indices."+str(repeat_idx)+"."+str(batch_size)+".pkl", "wb") as f:
            pickle.dump(testset.indices, f)

    return test_loader