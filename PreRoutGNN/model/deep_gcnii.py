import torch
import torch.nn.functional as F
import dgl
import dgl.function as fn
import functools
import pdb
from .mlp import MLP
from torch import nn
from utils import record_time

# {AllConv, DeepGCNII}: Simple and Deep Graph Convolutional Networks, arxiv 2007.02133 (GCNII)
class AllConv(torch.nn.Module):
    def __init__(self, in_nf, in_ef, out_nf, dropout, h=12):
        super().__init__()
        self.in_nf = in_nf
        self.in_ef = in_ef
        self.out_nf = out_nf
        self.h = h
        
        self.MLP_msg = MLP(in_nf*2+in_ef, 64, 64, 64, 1+4*h, dropout=dropout)
        self.MLP_reduce = MLP(in_nf+4*h, 64, 64, 64, out_nf, dropout=dropout)
        self.reset()
    
    def reset(self):
        return
        # for m in self.modules():
            # if isinstance(m, nn.Linear):
                # nn.init.xavier_normal_(m.weight)


    def edge_udf(self, edges):
        x = self.MLP_msg(torch.cat([edges.src['nf'], edges.dst['nf'], edges.data['ef']], dim=1))
        k, f1, f2, f3, f4 = torch.split(x, [1, self.h, self.h, self.h, self.h], dim=1)
        k = torch.sigmoid(k)
        return {
            "ef1": f1 * k,
            "ef2": f2 * k,
            "ef3": f3 * k,
            "ef4": f4 * k,
        }


    def forward(self, g, nf):
        with g.local_scope():
            g.ndata['nf'] = nf
            g.apply_edges(self.edge_udf)
            g.update_all(fn.copy_e('ef1', 'ef1'), fn.mean('ef1', 'nf1'))
            g.update_all(fn.copy_e('ef2', 'ef2'), fn.sum('ef2', 'nf2'))
            g.update_all(fn.copy_e('ef3', 'ef3'), fn.min('ef3', 'nf3'))
            g.update_all(fn.copy_e('ef4', 'ef4'), fn.max('ef4', 'nf4'))
            x = torch.cat([g.ndata['nf'], g.ndata['nf1'], g.ndata['nf2'], g.ndata['nf3'], g.ndata['nf4']], dim=1)
            x = self.MLP_reduce(x)
            return x


class DeepGCNII(torch.nn.Module):
    def __init__(self, in_nf, in_ef, out_nf, hidden_dim, dropout, n_layers):
        super().__init__()
        self.in_nf = in_nf
        self.in_ef = in_ef
        self.out_nf = out_nf
        self.hidden_dim = hidden_dim
        self.n_layers = n_layers
        
        self.layer0 = AllConv(in_nf, in_ef, hidden_dim, dropout)
        self.layers = [AllConv(in_nf+hidden_dim, in_ef, hidden_dim, dropout) for i in range(n_layers - 2)]
        self.layern = AllConv(hidden_dim, in_ef, out_nf, dropout)
        self.layers_store = torch.nn.Sequential(*self.layers)

    # @record_time
    def forward(self, g, nf):
        x = self.layer0(g, nf)
        for layer in self.layers:
            x = layer(g, torch.cat([x, nf], dim=1)) + x
        x = self.layern(g, x)
        return x
