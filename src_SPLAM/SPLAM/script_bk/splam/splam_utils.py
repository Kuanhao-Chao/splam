#! /usr/bin/python3
"""All the utility function for SPLAM!

    File name: splam_utils.py
    Author: Kuan-Hao Chao
    Email: kh.chao@cs.jhu.edu
    Date created: 12/20/2022
    Date last modified: 01/14/2023
    Python Version: 3.8
"""

import torch
from torch.optim.lr_scheduler import LambdaLR
from torch.optim import Optimizer
from torch.utils.data import Dataset, DataLoader, random_split

import numpy as np
import math
import random
import pickle
from sklearn.metrics import average_precision_score

JUNC_START = 200
JUNC_END = 600
IN_MAP = np.asarray([[0, 0, 0, 0],
                     [1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
# One-hot encoding of the inputs: 0 is for padding, and 1, 2, 3, 4 correspond
# to A, C, G, T respectively.
OUT_MAP = np.asarray([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1],
                      [0, 0, 0]])
# One-hot encoding of the outputs: 0 is for no splice, 1 is for acceptor,
# 2 is for donor and -1 is for padding.

           
def one_hot_encode_classifier(Xd):
    return IN_MAP[Xd.astype('int8')]

      
def create_datapoints(seq, strand):
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0').replace('K', '0').replace('R', '0').replace('Y', '0').replace('M', '0')
    jn_start = JUNC_START
    jn_end = JUNC_END

    #######################################
    # predicting splice / non-splice
    #######################################
    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = 0
    if strand == '+':
        Y0 = 1
    X = one_hot_encode_classifier(X0)
    return X, Y0


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


def get_accuracy(y_prob, y_true):
    assert y_true.ndim == 1 and y_true.size() == y_prob.size()
    y_prob = y_prob > 0.5
    return (y_true == y_prob).sum().item() / y_true.size(0)


def model_fn(DNAs, labels, model, criterion):
    """Forward a batch through the model."""
    outs = model(DNAs)
    outs = torch.flatten(outs)
    loss, accuracy = categorical_crossentropy_2d(labels, outs, criterion)
    return loss, accuracy, outs


def weighted_binary_cross_entropy(output, target, weights=None):    
    if weights is not None:
        assert len(weights) == 2
        loss = weights[1] * (target * torch.log(output+1e-10)) + \
               weights[0] * ((1 - target) * torch.log(1 - output+1e-10))
    else:
        loss = target * torch.log(output+1e-10) + (1 - target) * torch.log(1 - output+1e-10)
    return torch.neg(torch.mean(loss))


def categorical_crossentropy_2d(y_true, y_pred, criterion):
    weights = torch.FloatTensor([1.0, 1.0]) 
    return weighted_binary_cross_entropy(y_pred, y_true, weights), get_accuracy(y_pred, y_true)