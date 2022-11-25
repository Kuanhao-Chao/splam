import torch
import torch.nn as nn
import numpy as np
import sys
sys.path.append('./Conformer')
from CF import Conformer
import csv

class SpliceNN_Conformer(nn.Module):
    def __init__(self, config, d_model, state_num=1, dropout=0.1):
        super().__init__()
        # Project the dimension of features from that of input into d_model.
        self.prenet = nn.Linear(4, d_model)
        # TODO:
        #   Change Transformer to Conformer.
        #   https://arxiv.org/abs/2005.08100
        #   https://github.com/Masao-Someki/Conformer
        self.encoder_layer = Conformer(**config)
        # Project the the dimension of features from d_model into speaker nums.
        self.pred_layer = nn.Sequential(
        #   nn.Linear(d_model, d_model),
        #   nn.ReLU(),
          nn.Linear(d_model, 3),
          nn.Softmax(dim=2)
        )


    def forward(self, input):
        """
        args:
          input: (batch size, length, 40)
        return:
          out: (batch size, n_spks)
        """
        out = self.prenet(input)
        # out: (batch size, length, d_model)
        # out = out.permute(1, 0, 2)
        # out: (length, batch size, d_model)
        # For transformer: The encoder layer expect features in the shape of (length, batch size, d_model).
        # For conformer: The encoder layer expect features in the shape of (batch size, length, d_model).
        out = self.encoder_layer(out)
        # out: (length, batch size, d_model)
        # conformer out: (batch size, length, d_model)
        # out: (batch size, length, d_model)
        # mean pooling
        # print("out.size() before: ", out.size())
        # stats = out.mean(dim=1)
        # print("stats.size(): ", stats)
        out = self.pred_layer(out)
        # print("out.size(): ", out.size())
        # print("out: ", out)
        # out: (batch, n_spks)
        return out


def train_parse_args():
    
    """arguments"""

    config = {
        "data_dir": "./INPUTS",
        "batch_size": 32,
        "n_workers": 2,
        "valid_steps": 2000,
        "warmup_steps": 1000,
        "save_steps": 10000,
        "total_steps": 200000,
        "model_config":{
            "config1":{
                "d_model":80,
                "ff1_hsize": 2048,
                "ff1_dropout": 0.2,
                "n_head": 2,
                "mha_dropout": 0.2,
                "kernel_size": 3,
                "conv_dropout": 0.2,
                "ff2_hsize": 2048,
                "ff2_dropout": 0.2
            },
            "config2":{
                "d_model":80,
                "ff1_hsize": 2048,
                "ff1_dropout": 0.2,
                "n_head": 4,
                "mha_dropout": 0.2,
                "kernel_size": 3,
                "conv_dropout": 0.2,
                "ff2_hsize": 2048,
                "ff2_dropout": 0.2
            },
            "config3":{
                "d_model":80,
                "ff1_hsize": 4096,
                "ff1_dropout": 0.2,
                "n_head": 2,
                "mha_dropout": 0.2,
                "kernel_size": 3,
                "conv_dropout": 0.2,
                "ff2_hsize": 4096,
                "ff2_dropout": 0.2
            },
        },
        "model_path": {
            "config1":"./MODEL/Conformer/model96214.ckpt",
            "config2":"./MODEL/Conformer/model96142.ckpt",
            "config3":"./MODEL/Conformer/model96642.ckpt"},
    }

    return config