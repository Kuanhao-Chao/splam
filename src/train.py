import os
import json
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, random_split
from torch.nn.utils.rnn import pad_sequence
from torch.optim import Optimizer, AdamW
from torch.optim.lr_scheduler import LambdaLR
import random
import math
from tqdm import tqdm
from pathlib import Path
import sys
sys.path.append('./Conformer')
from CF import Conformer
import csv


class Classifier(nn.Module):
    def __init__(self, config, d_model=80, n_spks=600, dropout=0.1):
        super().__init__()
        # Project the dimension of features from that of input into d_model.
        self.prenet = nn.Linear(40, d_model)
        # TODO:
        #   Change Transformer to Conformer.
        #   https://arxiv.org/abs/2005.08100
        #   https://github.com/Masao-Someki/Conformer
        self.encoder_layer = Conformer(**config)
        #self.encoder_layer = nn.TransformerEncoderLayer(
        #  d_model=d_model, dim_feedforward=256, nhead=1
        #)
        # self.encoder = nn.TransformerEncoder(self.encoder_layer, num_layers=2)

        # Project the the dimension of features from d_model into speaker nums.
        self.pred_layer = nn.Sequential(
          #nn.Linear(d_model, d_model),
          #nn.ReLU(),
          nn.Linear(d_model, n_spks),
        )


    def forward(self, mels):
        """
        args:
          mels: (batch size, length, 40)
        return:
          out: (batch size, n_spks)
        """
        out = self.prenet(mels)
        # out: (batch size, length, d_model)
        #out = out.permute(1, 0, 2)
        # out: (length, batch size, d_model)
        # For transformer: The encoder layer expect features in the shape of (length, batch size, d_model).
        # For conformer: The encoder layer expect features in the shape of (batch size, length, d_model).
        out = self.encoder_layer(out)
        # out: (length, batch size, d_model)
        # conformer out: (batch size, length, d_model)
        #out = out.transpose(0, 1)
        # out: (batch size, length, d_model)
        # mean pooling
        stats = out.mean(dim=1)
        # stats: (batch size, d_model)
        out = self.pred_layer(stats)
        # out: (batch, n_spks)
        return out



