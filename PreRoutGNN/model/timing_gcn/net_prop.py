from model import MLP
import torch
import torch.nn.functional as F
import dgl
import dgl.function as fn
import functools
from config import Config


class NetConv(torch.nn.Module):
    """return new node embedding"""
    def __init__(self, in_nf, in_ef, out_nf, dropout, h1=32, h2=32):
        super().__init__()
        self.in_nf = in_nf
        self.in_ef = in_ef
        self.out_nf = out_nf
        self.h1 = h1
        self.h2 = h2
        self.hidden_dim = Config.hidden_dim_netprop
        
        self.MLP_msg_i2o = MLP(self.in_nf * 2 + self.in_ef, self.hidden_dim, self.hidden_dim, self.hidden_dim, 1 + self.h1 + self.h2, dropout=dropout)
        self.MLP_reduce_o = MLP(self.in_nf + self.h1 + self.h2, self.hidden_dim, self.hidden_dim, self.hidden_dim, self.out_nf, dropout=dropout)
        self.MLP_msg_o2i = MLP(self.in_nf * 2 + self.in_ef, self.hidden_dim, self.hidden_dim, self.hidden_dim, self.hidden_dim, self.out_nf, dropout=dropout)

    def edge_msg_i(self, edges):
        '''use edge src node embedding, edge dst node embedding, edge embedding to generate message'''
        x = torch.cat([edges.src['nf'], edges.dst['nf'], edges.data['ef']], dim=1)
        x = self.MLP_msg_o2i(x)
        return {'efi': x}

    def edge_msg_o(self, edges):
        '''use edge src node embedding, edge dst node embedding, edge embedding to generate 2 messages'''
        x = torch.cat([edges.src['nf'], edges.dst['nf'], edges.data['ef']], dim=1)
        x = self.MLP_msg_i2o(x)
        k, f1, f2 = torch.split(x, [1, self.h1, self.h2], dim=1)
        k = torch.sigmoid(k)
        return {'efo1': f1 * k, 'efo2': f2 * k}

    def node_reduce_o(self, nodes):
        x = torch.cat([nodes.data['nf'], nodes.data['nfo1'], nodes.data['nfo2']], dim=1)
        x = self.MLP_reduce_o(x)
        return {'new_nf': x}
        
    def forward(self, g, ts, nf):
        with g.local_scope():
            g.ndata['nf'] = nf
            g.update_all(self.edge_msg_i, fn.sum('efi', 'new_nf'), etype='net_out')
            g.apply_edges(self.edge_msg_o, etype='net_in')
            g.update_all(fn.copy_e('efo1', 'efo1'), fn.sum('efo1', 'nfo1'), etype='net_in')
            g.update_all(fn.copy_e('efo2', 'efo2'), fn.max('efo2', 'nfo2'), etype='net_in')
            g.apply_nodes(self.node_reduce_o, ts['output_nodes'])
            
            return g.ndata['new_nf']