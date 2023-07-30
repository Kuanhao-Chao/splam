from torch.nn import Module, BatchNorm1d, LazyBatchNorm1d, ReLU, LeakyReLU, Conv1d, LazyConv1d, ModuleList, Softmax, Sigmoid, Flatten, Dropout2d, Linear
import numpy as np

CARDINALITY_ITEM = 16

class ResidualUnit(Module):
    def __init__(self, l, w, ar, bot_mul=1):
        super().__init__()
        bot_channels = int(round(l * bot_mul))
        self.batchnorm1 = BatchNorm1d(l)
        self.relu = LeakyReLU(0.1)
        self.batchnorm2 = BatchNorm1d(l)
        self.C = bot_channels//CARDINALITY_ITEM
        self.conv1 = Conv1d(l, l, w, dilation=ar, padding=(w-1)*ar//2, groups=self.C)
        self.conv2 = Conv1d(l, l, w, dilation=ar, padding=(w-1)*ar//2, groups=self.C)

    def forward(self, x, y):
        x1 = self.relu(self.batchnorm1(self.conv1(x)))
        x2 = self.relu(self.batchnorm2(self.conv2(x1)))
        # x1 = self.relu(self.batchnorm1(self.conv1(x)))
        # x2 = self.relu(self.batchnorm1(self.conv1(x1)))

        # print("x : ", x.size())
        # print("x1: ", x1.size())
        # print("x2: ", x2.size())
        return x + x2, y


class Skip(Module):
    def __init__(self, l):
        super().__init__()
        self.conv = Conv1d(l, l, 1)

    def forward(self, x, y):
        return x, self.conv(x) + y


class SPLAM(Module):
    def __init__(self, L=64, W=np.array([11]*8+[21]*4+[41]*4), AR=np.array([1]*4+[4]*4+[10]*4+[25]*4)):
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
        self.last_cov = Conv1d(L, 3, 1)
        self.softmax = Softmax(dim=1)
        # self.flatten = Flatten()
        # self.drop_out = Dropout2d(0.2)
        # self.fc = Linear(2400, 1)
        # self.softmax = Softmax(dim=1)
        # self.sigmoid = Sigmoid()

    def forward(self, x):
        x, skip = self.skip1(self.conv1(x), 0)
        for m in self.residual_blocks:
            x, skip = m(x, skip)
        #######################################
        # predicting pb for every bp
        #######################################
        return self.softmax(self.last_cov(skip))

        # #######################################
        # # predicting splice / non-splice
        # #######################################
        # output = self.sigmoid(self.fc(self.flatten(self.last_cov(skip))))
        # return output