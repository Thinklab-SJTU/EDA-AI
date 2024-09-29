from model import MLP
from config import Config
from einops import rearrange
from torch import nn

import torch
import torch.nn.functional as F
import dgl
import dgl.function as fn
import functools
 


class SignalPropAttn5(torch.nn.Module):
    """
    add residual connection in edge_msg_net_prop
    based on ResidualSignalProp4
    num_keys = k1 * k2, which is 7*7 == 49
    only one query
    """
    def __init__(self, in_nf, in_cell_num_luts, in_cell_lut_sz, out_nf, out_cef, dropout, h1=32, h2=32, lut_dup=4):
        super().__init__()
        self.in_nf = in_nf
        self.in_cell_num_luts = in_cell_num_luts
        self.in_cell_lut_sz = in_cell_lut_sz
        self.out_nf = out_nf
        self.out_cef = out_cef
        self.h1 = h1
        self.h2 = h2
        self.lut_dup = lut_dup
        self.num_heads = self.lut_dup
        self.hidden_dim_cellprop = Config.hidden_dim_cellprop
        
        self.qdim = 32
        self.kdim = 1
        self.vdim = 1
        self.msa = nn.MultiheadAttention(embed_dim=self.qdim, kdim=self.kdim, vdim=self.vdim, batch_first=True, dropout=dropout, num_heads=self.num_heads)
        
        # generate net prop edge message
        self.MLP_netprop = MLP(self.out_nf + 2 * self.in_nf, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.out_nf, dropout=dropout)
        
        # generate cell prop edge message
        self.MLP_lut_query = MLP(self.out_nf + 2 * self.in_nf, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.in_cell_num_luts * self.qdim, dropout=dropout)
        self.MLP_cellarc_msg = MLP(self.out_nf + 2 * self.in_nf + self.qdim*self.in_cell_num_luts, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.hidden_dim_cellprop, 1 + self.h1 + self.h2 + self.out_cef, dropout=dropout)
        
        # reduce messages during cell prop
        self.MLP_cellreduce = MLP(self.in_nf + self.h1 + self.h2, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.hidden_dim_cellprop, self.out_nf, dropout=dropout)


    def node_reduce_primary_input(self, nodes):
        return {'new_nf': nodes.data['n_atslew']}


    def edge_msg_net_prop(self, edges, groundtruth=False):
        if groundtruth:
            last_nf = edges.src['n_atslew']
        else:
            last_nf = edges.src['new_nf']
        
        x = torch.cat([last_nf, edges.src['nf'], edges.dst['nf']], dim=1)
        x = self.MLP_netprop(x)
        return {'efn': x + last_nf}


    def edge_msg_cell_prop(self, edges, groundtruth=False):
        if groundtruth:
            last_nf = edges.src['n_atslew']
        else:
            last_nf = edges.src['new_nf']
        
        # [EL, 1, qdim]
        q = torch.cat([last_nf, edges.src['nf'], edges.dst['nf']], dim=1)
        q = self.MLP_lut_query(q)
        q = q.reshape(-1, 1, self.qdim)

        # [EL, 49, kdim]
        axis_len = self.in_cell_num_luts * (1 + 2 * self.in_cell_lut_sz)
        axis = edges.data['ef'][:, :axis_len]
        axis = rearrange(axis, 'E (L D) -> E L D', L=self.in_cell_num_luts, D=1+2*self.in_cell_lut_sz)
        key1 = axis[:, :, 1:1+self.in_cell_lut_sz]
        key2 = axis[:, :, 1+self.in_cell_lut_sz:1+self.in_cell_lut_sz+self.in_cell_lut_sz]
        k = torch.einsum('e l i, e l j -> e l i j', key1, key2)
        k = rearrange(k, 'e l i j -> (e l) (i j) 1')
        
        # [EL, 49, vdim]
        tables_len = self.in_cell_num_luts * self.in_cell_lut_sz ** 2
        tables = edges.data['ef'][:, axis_len:axis_len + tables_len]
        v = rearrange(tables, 'E (L D) -> (E L) D 1', L=self.in_cell_num_luts, D=self.in_cell_lut_sz**2)
        
        # [EL, 1, qdim]
        output = self.msa(q,k,v)[0]
        r = output.reshape(-1, self.qdim*self.in_cell_num_luts)
        x = torch.cat([last_nf, edges.src['nf'], edges.dst['nf'], r], dim=1)
        
        x = self.MLP_cellarc_msg(x)
        k, f1, f2, cef = torch.split(x, [1, self.h1, self.h2, self.out_cef], dim=1)
        k = torch.sigmoid(k)
        return {'efc1': f1 * k, 'efc2': f2 * k, 'efce': cef}


    def node_reduce_cell_prop(self, nodes):
        x = torch.cat([nodes.data['nf'], nodes.data['nfc1'], nodes.data['nfc2']], dim=1)
        x = self.MLP_cellreduce(x)
        return {'new_nf': x}

        
    def forward(self, g, ts, nf, groundtruth=False):
        assert len(ts['topo']) % 2 == 0, 'The number of logic levels must be even (net, cell, net)'
        
        with g.local_scope():
            g.ndata['nf'] = nf
            g.ndata['new_nf'] = torch.zeros(g.num_nodes(), self.out_nf, device=g.device, dtype=nf.dtype)
            
            g.apply_nodes(self.node_reduce_primary_input, ts['pi_nodes'])

            def prop_net(nodes, groundtruth):
                g.pull(nodes, functools.partial(self.edge_msg_net_prop, groundtruth=groundtruth), fn.sum('efn', 'new_nf'), etype='net_out')

            def prop_cell(nodes, groundtruth):

                es = g.in_edges(nodes, etype='cell_out')
                g.apply_edges(functools.partial(self.edge_msg_cell_prop, groundtruth=groundtruth), es, etype='cell_out')
                g.send_and_recv(es, fn.copy_e('efc1', 'efc1'), fn.sum('efc1', 'nfc1'), etype='cell_out')
                g.send_and_recv(es, fn.copy_e('efc2', 'efc2'), fn.max('efc2', 'nfc2'), etype='cell_out')
                g.apply_nodes(self.node_reduce_cell_prop, nodes)
            
            if groundtruth:
                prop_net(ts['input_nodes'], groundtruth)
                prop_cell(ts['output_nodes_nonpi'], groundtruth)

            else:
                for i in range(1, len(ts['topo'])):
                    if i % 2 == 1:
                        prop_net(ts['topo'][i], groundtruth)
                    else:
                        prop_cell(ts['topo'][i], groundtruth)
            
            return g.ndata['new_nf'], g.edges['cell_out'].data['efce']
        