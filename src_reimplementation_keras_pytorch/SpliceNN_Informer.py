import os
import json
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, random_split
from torch.nn.utils.rnn import pad_sequence
from torch.optim import Optimizer, AdamW
from torch.optim.lr_scheduler import LambdaLR
import math
from tqdm import tqdm
from pathlib import Path
import sys
from models.model import Informer, InformerStack
# from longformer.longformer import Longformer
# from transformers import LongformerModel, LongformerConfig
# from memory_transformer_xl import MemoryTransformerXL
# from CF import Conformer
import csv

class Informer_SpliceAI(nn.Module):
    def __init__(self, config, seq_len, d_model=10, dropout=0.1):
        super().__init__()
        # Project the dimension of features from that of input into d_model.
        self.prenet = nn.Linear(4, d_model)
        #   Change Transformer to Conformer.
        #   https://arxiv.org/abs/2005.08100
        #   https://github.com/Masao-Someki/Conformer
        self.informer = Informer(**config)

        # Project the the dimension of features from d_model into speaker nums.
        self.pred_layer = nn.Sequential(
          # nn.Linear(d_model, d_model),
          #nn.ReLU(),
          nn.Linear(d_model, 3),
          nn.Softmax(dim=2)
        )

    def forward(self, DNAs):
        """
        args:
          DNAs: (batch size, length, 4)
        return:
          out: (batch size, length, 3)
        """
        # print("DNAs.shape: ", DNAs.size())
        out = self.prenet(DNAs)
        # out: (batch size, length, d_model)
        #out = out.permute(1, 0, 2)
        # out: (length, batch size, d_model)
        # For transformer: The encoder layer expect features in the shape of (length, batch size, d_model).
        # For conformer: The encoder layer expect features in the shape of (batch size, length, d_model).
        out = self.informer(out)
        # out = self.encoder_layer(out)
        # out: (length, batch size, d_model)
        # conformer out: (batch size, length, d_model)
        #out = out.transpose(0, 1)
        # out: (batch size, length, d_model)
        # mean pooling
        # stats = out.mean(dim=1)
        # print("stats.size(): ", stats.size())
        # stats: (batch size, d_model)
        out = self.pred_layer(out)
        # print("out: ", out)
        # out: (batch, n_spks)
        return out


def get_cosine_schedule_with_warmup(
      optimizer: Optimizer,
      num_warmup_steps: int,
      num_training_steps: int,
      num_cycles: float = 0.5,
      last_epoch: int = -1,
    ):
    """
    Create a schedule with a learning rate that decreases following the values of the cosine function between the
    initial lr set in the optimizer to 0, after a warmup period during which it increases linearly between 0 and the
    initial lr set in the optimizer.

    Args:
    optimizer (:class:`~torch.optim.Optimizer`):
      The optimizer for which to schedule the learning rate.
    num_warmup_steps (:obj:`int`):
      The number of steps for the warmup phase.
    num_training_steps (:obj:`int`):
      The total number of training steps.
    num_cycles (:obj:`float`, `optional`, defaults to 0.5):
      The number of waves in the cosine schedule (the defaults is to just decrease from the max value to 0
      following a half-cosine).
    last_epoch (:obj:`int`, `optional`, defaults to -1):
      The index of the last epoch when resuming training.

    Return:
    :obj:`torch.optim.lr_scheduler.LambdaLR` with the appropriate schedule.
    """

    def lr_lambda(current_step):
        # Warmup
        if current_step < num_warmup_steps:
            return float(current_step) / float(max(1, num_warmup_steps))
        # decadence
        progress = float(current_step - num_warmup_steps) / float(
          max(1, num_training_steps - num_warmup_steps)
        )
        return max(
          0.0, 0.5 * (1.0 + math.cos(math.pi * float(num_cycles) * 2.0 * progress))
        )

    return LambdaLR(optimizer, lr_lambda, last_epoch)


def categorical_crossentropy_2d(y_true, y_pred):
    # prod = output[:,0]*target
    # return -prod[prod<0].sum()

    return - torch.mean(1*y_true[:, :, 0]*torch.log(y_pred[:, :, 0]+1e-10)
                        + 1000*y_true[:, :, 1]*torch.log(y_pred[:, :, 1]+1e-10)
                        + 1000*y_true[:, :, 2]*torch.log(y_pred[:, :, 2]+1e-10))


def model_fn(DNAs, labels, model, criterion, device):
    """Forward a batch through the model."""

    # DNAs, labels = batch # size [batch size, length, 4] , [batch size, length]

    # print("\nDNAs: ", DNAs.size())
    # print("labels: ", labels.size())
    # DNAs = DNAs.to(torch.float32).to(device)
    # labels = labels.to(torch.float32).to(device)

    # labels = torch.stack(list(labels), dim=0)
    # print("\nmels: ", mels.size())
    # print("labels hah: ", labels.size())
    # print("labels len: ", len(labels))
    # print("labels: ", labels[0].size())

    outs = model(DNAs)
    # print("outs.size(): ", outs.size())
    # outs: (batch, n_spks)

    loss = categorical_crossentropy_2d(outs, labels)
    # criterion(outs, labels)

    # # Get the speaker id with highest probability.
    # preds = (outs>0.5).float()
    # # print("preds: ", preds)
    # # print("preds.size(): ", preds.size())
    # # Compute accuracy.
    # accuracy = torch.mean((preds == labels).float())

    return loss, outs




# def valid(dataloader, model, criterion, device):
#     """Validate on validation set."""

#     model.eval()
#     running_loss = 0.0
#     running_accuracy = 0.0
#     pbar = tqdm(total=len(dataloader.dataset), ncols=0, desc="Valid", unit=" uttr")

#     for i, batch in enumerate(dataloader):
#         with torch.no_grad():
#             loss, accuracy = model_fn(batch, model, criterion, device)
#             running_loss += loss.item()
#             running_accuracy += accuracy.item()

#         pbar.update(dataloader.batch_size)
#         pbar.set_postfix(
#           loss=f"{running_loss / (i+1):.2f}",
#           accuracy=f"{running_accuracy / (i+1):.2f}",
#         )

#     pbar.close()
#     model.train()

#     return running_accuracy / len(dataloader)





def train_parse_args():

    """arguments"""

    config = {
        "data_dir": "../Dataset",
        "batch_size": 32,
        "n_workers": 8,
        "valid_steps": 2000,
        "warmup_steps": 1000,
        "save_steps": 10000,
        "total_steps": 200000,
        "model_config":{
            "config1":{
                "enc_in": 10, 
                "dec_in": 10, 
                "c_out": 3, 
                "seq_len": 5000, 
                "label_len": 48, 
                "out_len": 3, 
                "factor": 5, 
                "d_model": 512,
                "n_heads": 8, 
                "e_layers": 3,
                "d_layers": 2,
                "d_ff": 512, 
                "dropout": 0.0, 
                "attn": 'prob',
                "embed": 'fixed',
                "freq": 'h',
                "activation": 'gelu', "output_attention": False,
                "distil": True,
                "mix": True,
            },
        },
        "model_path": {
            "config1":"./model96214.ckpt",
            "config2":"./model96142.ckpt",
            "config3":"./model96642.ckpt"},
    }

    return config

# fix random seed
def same_seeds(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True



class InferenceDataset(Dataset):
    def __init__(self, data_dir):
        testdata_path = Path(data_dir) / "testdata.json"
        metadata = json.load(testdata_path.open())
        self.data_dir = data_dir
        self.data = metadata["utterances"]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        utterance = self.data[index]
        feat_path = utterance["feature_path"]
        mel = torch.load(os.path.join(self.data_dir, feat_path))

        return feat_path, mel


def inference_collate_batch(batch):
    """Collate a batch of data."""
    feat_paths, mels = zip(*batch)

    return feat_paths, torch.stack(mels)

if __name__ == "__main__":

    train_main(**train_parse_args())
    test_main(**test_parse_args())



# References
# This code is modified from TA's sample code in NTU machine learning course
# Conformer source I used: https://github.com/Masao-Someki/Conformer
