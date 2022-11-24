from torch.nn import Module, BatchNorm1d, ReLU, Conv1d, ModuleList, Softmax
import numpy as np


CL_max = 10000


class ResidualUnit(Module):
    def __init__(self, l, w, ar):
        super().__init__()
        self.batchnorm1 = BatchNorm1d(l)
        self.relu = ReLU()
        self.batchnorm2 = BatchNorm1d(l)
        self.conv1 = Conv1d(l, l, w, dilation=ar, padding=(w-1)*ar//2)
        self.conv2 = Conv1d(l, l, w, dilation=ar, padding=(w-1)*ar//2)

    def forward(self, x, y):
        x1 = self.conv1(self.relu(self.batchnorm1(x)))
        x2 = self.conv2(self.relu(self.batchnorm2(x1)))
        return x + x2, y


class Skip(Module):
    def __init__(self, l):
        super().__init__()
        self.conv = Conv1d(l, l, 1)

    def forward(self, x, y):
        return x, self.conv(x) + y


class SpliceNN(Module):
    def __init__(self, L=32, W=np.array([11]*8+[21]*4+[41]*4), AR=np.array([1]*4+[4]*4+[10]*4+[25]*4)):
        super().__init__()
        self.CL = 2 * (AR * (W - 1)).sum()  # context length
        self.conv1 = Conv1d(4, L, 1)
        self.skip1 = Skip(L)
        self.residual_blocks = ModuleList()
        for i, (w, r) in enumerate(zip(W, AR)):
            self.residual_blocks.append(ResidualUnit(L, w, r))
            if (i+1) % 4 == 0:
                self.residual_blocks.append(Skip(L))
        if (len(W)+1) % 4 != 0:
            self.residual_blocks.append(Skip(L))
        self.last = Conv1d(L, 3, 1)
        self.softmax = Softmax(dim=1)

    def forward(self, x):
        x, skip = self.skip1(self.conv1(x), 0)
        for m in self.residual_blocks:
            x, skip = m(x, skip)
        return self.softmax(self.last(skip)[..., CL_max//2:-CL_max//2])
