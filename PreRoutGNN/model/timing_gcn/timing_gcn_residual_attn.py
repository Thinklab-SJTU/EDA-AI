from model import MLP
from .net_prop import NetConv
from .cell_prop_attn import *
from torch import nn
from config import Config
import torch
import torch.nn.functional as F
import dgl
import dgl.function as fn
import functools
from utils import record_time



class PreRoutGNN(torch.nn.Module):
    """
    Add residual on net conv.
    Selecting 4 components as netdelay -> use a linear to get 4 components as netdelay.
    """
    def __init__(self, in_nf, out_nf, dropout, *args, **kwds):
        super().__init__()
        self.in_nf = in_nf
        self.ef = 2 if Config.num_pin_location_freq_compents <= 0 else 2 + 2*2*Config.num_pin_location_freq_compents
        self.hidden_dim = Config.hidden_dim_gcn
        self.out_nf = out_nf

        self.dim_up = nn.Sequential(nn.Linear(self.in_nf, self.hidden_dim), nn.LeakyReLU(0.2))
        self.to_netdelay = MLP(self.hidden_dim, self.hidden_dim, 4, dropout=dropout)

        self.nc1 = NetConv(self.hidden_dim, self.ef, self.hidden_dim, dropout)
        self.nc2 = NetConv(self.hidden_dim, self.ef, self.hidden_dim, dropout)
        self.nc3 = NetConv(self.hidden_dim, self.ef, self.hidden_dim, dropout)  
        self.nc4 = NetConv(self.hidden_dim, self.ef, self.hidden_dim, dropout)  

        self.prop = SignalPropAttn5(self.in_nf + self.hidden_dim, 8, 7, self.out_nf, 4, dropout)

    # @record_time
    def forward(self, g, ts, groundtruth=False):
        nf0 = g.ndata['nf']
        x = self.dim_up(nf0)
        x = self.nc1(g, ts, x) + x
        x = self.nc2(g, ts, x) + x
        x = self.nc3(g, ts, x) + x
        x = self.nc4(g, ts, x) + x

        # net delay is pin/node label
        net_delays = self.to_netdelay(x)

        nf1 = torch.cat([nf0, x], dim=1)
        atslew, cell_delays = self.prop(g, ts, nf1, groundtruth=groundtruth)
        return net_delays, cell_delays, atslew