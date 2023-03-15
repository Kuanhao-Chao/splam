import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
import random

class myDataset(Dataset):
    def __init__(self, X, Y, segment_len=5000):
        self.segment_len = segment_len
        self.X = X
        self.Y = Y
        # print("len(Y): ", len(self.X))
 
    def __len__(self):
        return len(self.X)
 
    def __getitem__(self, idx):
        # print("X: ", self.X.shape)
        # print("Y: ", self.Y[0].shape)
        x = self.X[idx]
        y = self.Y[0][idx]
        x = torch.FloatTensor(x).long()
        y = torch.FloatTensor(y).long()
        # print("x: ", x.size())
        # print("y: ", y.size())
        return x, y

