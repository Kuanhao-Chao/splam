import torch
from torch.nn import Module, BatchNorm1d, LeakyReLU, Conv1d, ModuleList, Softmax, Sigmoid, Flatten, Dropout2d, Linear
from torch.optim.lr_scheduler import LambdaLR
from torch.optim import Optimizer
from torch.utils.data import Dataset, DataLoader

import math
import random
import numpy as np
import warnings
import os
import re
from progress.bar import Bar

warnings.filterwarnings('ignore')

CARDINALITY_ITEM = 16
SEQ_LEN = 800
JUNC_START = 200
JUNC_END = 600
IN_MAP = np.asarray([[0, 0, 0, 0],
                     [1, 0, 0, 0],
                     [0, 1, 0, 0],
                     [0, 0, 1, 0],
                     [0, 0, 0, 1]])
OUT_MAP = np.asarray([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1],
                      [0, 0, 0]])

def get_donor_acceptor_scores(D_YL, A_YL, D_YP, A_YP):
    return D_YL[:, 200], D_YP[:, 200], A_YL[:, 600], A_YP[:, 600]

def one_hot_encode(Xd, Yd):
    return IN_MAP[Xd.astype('int8')], [OUT_MAP[Yd[t].astype('int8')] for t in range(1)]

def create_datapoints(seq, strand):
    # seq = 'N'*(CL_MAX//2) + seq + 'N'*(CL_MAX//2)
    seq = seq.upper().replace('A', '1').replace('C', '2').replace('G', '3').replace('T', '4')
    pattern = r'[^1234]'
    # Replace non-ACGT characters with 0
    seq = re.sub(pattern, '0', seq)
    jn_start = JUNC_START
    jn_end = JUNC_END

    #######################################
    # predicting pb for every bp
    #######################################
    X0 = np.asarray(list(map(int, list(seq))))
    Y0 = [np.zeros(SEQ_LEN) for t in range(1)]
    if strand == '+':
        for t in range(1):
            Y0[t][jn_start] = 2
            Y0[t][jn_end] = 1
    X, Y = one_hot_encode(X0, Y0)
    return X, Y

def get_cosine_schedule_with_warmup(
      optimizer: Optimizer,
      num_warmup_steps: int,
      num_training_steps: int,
      num_cycles: float = 0.5,
      last_epoch: int = -1,
    ):
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
    outs = model(DNAs)
    loss = categorical_crossentropy_2d(labels, outs, criterion)
    return loss, outs

def categorical_crossentropy_2d(y_true, y_pred, criterion):
    SEQ_WEIGHT = 5
    gamma = 2
    return - torch.mean(y_true[:, 0, :] * torch.mul( torch.pow( torch.sub(1, y_pred[:, 0, :]), gamma ), torch.log(y_pred[:, 0, :]+1e-10) )
                        + SEQ_WEIGHT * y_true[:, 1, :] * torch.mul( torch.pow( torch.sub(1, y_pred[:, 1, :]), gamma ), torch.log(y_pred[:, 1, :]+1e-10) )
                        + SEQ_WEIGHT * y_true[:, 2, :] * torch.mul( torch.pow( torch.sub(1, y_pred[:, 2, :]), gamma ), torch.log(y_pred[:, 2, :]+1e-10) ))

def split_seq_name(seq):
    return seq[1:]

class myDataset(Dataset):
    def __init__(self, type, of, shuffle, segment_len=800):
        self.segment_len = segment_len
        self.data = []
        self.indices = []
        pidx = 0
        with open(of, 'r') as f:
            lines = f.read().splitlines()
            seq_name = ''
            seq = ''
            for line in lines:
                if pidx % 2 == 0:
                    seq_name = split_seq_name(line)
                elif pidx % 2 == 1:
                    seq = line
                    if seq[0] == '>':
                        seq_name = line
                        continue
                    X, Y = create_datapoints(seq, '+')
                    X = torch.Tensor(np.array(X))
                    Y = torch.Tensor(np.array(Y)[0])
                    if X.size()[0] != 800:
                        print('seq_name: ', seq_name)
                        print(X.size())
                        print(Y.size())
                    self.data.append([X, Y, seq_name])
                pidx += 1
                if pidx %100000 == 0:
                    print('\t', pidx//2, ' junctions loaded.')

        index_shuf = list(range(len(self.data)))
        if shuffle:
            random.shuffle(index_shuf)
        list_shuf = [self.data[i] for i in index_shuf]
        self.data = list_shuf
        self.indices = index_shuf
        print('\t', pidx//2, ' junctions loaded.')

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        feature = self.data[index][0]
        label = self.data[index][1]
        seq_name = self.data[index][2]
        feature = torch.flatten(feature, start_dim=1)
        return feature, label, seq_name

def get_dataloader(batch_size, n_workers, output_file, shuffle, repeat_idx):
    testset = myDataset('test', output_file, shuffle, SEQ_LEN)
    test_loader = DataLoader(
        testset,
        batch_size = batch_size,
        drop_last = False,
        pin_memory = True,
    )
    return test_loader

def splam_prediction(junction_fasta, out_score_f, model_path, batch_size, device_str):
    BATCH_SIZE = int(batch_size)
    N_WORKERS = None
    if device_str == "NONE":
        device_str = 'cpu'
        if torch.cuda.is_available(): 
            device_str = 'cuda'
        elif torch.backends.mps.is_available():
            device_str = 'mps'
    device = torch.device(device_str)

    print(f'[Info] Running model in "'+ device_str+'" mode')
    print(f'[Info] Loading model ... (' + model_path + ')', flush = True)
    # model = torch.jit.load(model_path)
    print("model = torch.load(model_path)!!")
    model = torch.load(model_path)
    model = model.to(device)

    print(f'[Info] Done loading model', flush = True)
    print(f'[Info] Loading data ...', flush = True)
    test_loader = get_dataloader(BATCH_SIZE, N_WORKERS, junction_fasta, False, str(0))
    print(f'[Info] Done loading data', flush = True)
    
    criterion = torch.nn.BCELoss()
    fw_junc_scores = open(out_score_f, 'w')

    model.eval()
    junc_counter = 0
    pbar = Bar('[Info] SPLAM! ', max=len(test_loader))
    with torch.no_grad():
        for batch_idx, data in enumerate(test_loader):
            # DNAs:  torch.Size([40, 800, 4])
            # labels:  torch.Size([40, 1, 800, 3])
            DNAs, labels, seqname = data
            DNAs = DNAs.to(torch.float32).to(device)
            labels = labels.to(torch.float32).to(device)

            DNAs = torch.permute(DNAs, (0, 2, 1))
            labels = torch.permute(labels, (0, 2, 1))
            loss, yp = model_fn(DNAs, labels, model, criterion)
            is_expr = (labels.sum(axis=(1,2)) >= 1)
            A_YL = labels[is_expr, 1, :].to('cpu').detach().numpy()
            A_YP = yp[is_expr, 1, :].to('cpu').detach().numpy()
            D_YL = labels[is_expr, 2, :].to('cpu').detach().numpy()
            D_YP = yp[is_expr, 2, :].to('cpu').detach().numpy()

            donor_labels, donor_scores, acceptor_labels, acceptor_scores = get_donor_acceptor_scores(D_YL, A_YL, D_YP, A_YP)
            for idx in range(len(yp)):

                eles = seqname[idx].split(';')
                if len(eles) == 7:
                    chr, start, end, strand, name, aln_num, trans = eles
                    if strand == '+':
                        fw_junc_scores.write(f'{chr}\t{str(start)}\t{str(end)}\t{name}\t{str(aln_num)}\t{strand}\t{str(donor_scores[idx])}\t{str(acceptor_scores[idx])}\t{trans}\n')
                    elif strand == '-':
                        fw_junc_scores.write(f'{chr}\t{str(end)}\t{str(start)}\t{name}\t{str(aln_num)}\t{strand}\t{str(donor_scores[idx])}\t{str(acceptor_scores[idx])}\t{trans}\n')

                else:
                    chr, start, end, strand, name, aln_num = eles
                    if strand == '+':
                        fw_junc_scores.write(f'{chr}\t{str(start)}\t{str(end)}\t{name}\t{str(aln_num)}\t{strand}\t{str(donor_scores[idx])}\t{str(acceptor_scores[idx])}\n')
                    elif strand == '-':
                        fw_junc_scores.write(f'{chr}\t{str(end)}\t{str(start)}\t{name}\t{str(aln_num)}\t{strand}\t{str(donor_scores[idx])}\t{str(acceptor_scores[idx])}\n')
                
                junc_counter += 1            
            # increment the progress bar
            pbar.next()

    pbar.finish()
    fw_junc_scores.close()
    return out_score_f