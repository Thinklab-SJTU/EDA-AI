import torch
from torch import nn

class L1_norm(nn.Module):
    def __init__(self):
        super(L1_norm, self).__init__()

    def forward(self,x):
        return torch.abs(x).sum().cuda()