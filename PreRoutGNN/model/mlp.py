import torch
from torch import nn

class MLP(nn.Module):
    def __init__(self, *sizes, dropout:float, batchnorm=False):
        super().__init__()
        fcs = []
        for i in range(1, len(sizes)):
            fcs.append(nn.Linear(sizes[i - 1], sizes[i]))
            if i < len(sizes) - 1:
                fcs.append(nn.LeakyReLU(negative_slope=0.2))

                if dropout > 0.0: fcs.append(nn.Dropout(p=dropout))
                if batchnorm: fcs.append(nn.BatchNorm1d(sizes[i]))

        self.layers = nn.Sequential(*fcs)

    def forward(self, x):
        return self.layers(x)
