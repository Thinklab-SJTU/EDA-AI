from torch import nn
from .mlp import MLP
import torch
import torch.nn.functional as F
import dgl
import dgl.function as fn


# {AllConv, DeepGCNII}: Simple and Deep Graph Convolutional Networks, arxiv 2007.02133 (GCNII)
class AllConv(torch.nn.Module):
    def __init__(self, in_nf, in_ef, out_nf, dropout, h=12):
        super().__init__()
        self.in_nf = in_nf
        self.in_ef = in_ef
        self.out_nf = out_nf
        self.h = h
        
        self.MLP_msg = MLP(in_nf*2+in_ef, 64, 64, 64, 1+2*h, dropout=dropout)
        self.MLP_reduce = MLP(in_nf+2*h, 64, 64, 64, out_nf, dropout=dropout)
    

    def edge_udf(self, edges):
        x = self.MLP_msg(torch.cat([edges.src['nf'], edges.dst['nf'], edges.data['ef']], dim=1))
        k, f1, f2= torch.split(x, [1, self.h, self.h], dim=1)
        k = torch.sigmoid(k)
        return {
            "ef1": f1 * k,
            "ef2": f2 * k,
        }


    def forward(self, g, nf):
        with g.local_scope():
            g.ndata['nf'] = nf
            g.apply_edges(self.edge_udf)
            g.update_all(fn.copy_e('ef1', 'ef1'), fn.max('ef1', 'nf1'))
            g.update_all(fn.copy_e('ef2', 'ef2'), fn.sum('ef2', 'nf2'))
            x = torch.cat([g.ndata['nf'], g.ndata['nf1'], g.ndata['nf2']], dim=1)
            x = self.MLP_reduce(x)
            return x


class GCN(torch.nn.Module):
    def __init__(
        self,
        input_node_dim,
        input_edge_dim,
        hidden_dim,
        output_node_dim,
        num_layers,
        dropout,
    ):
        super().__init__()
        self.input_node_dim = input_node_dim
        self.input_edge_dim = input_edge_dim
        self.output_node_dim = output_node_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.dropout = dropout
        
        self.layer0 = AllConv(input_node_dim, input_edge_dim, hidden_dim, dropout)
        self.layers = [AllConv(input_node_dim+hidden_dim, input_edge_dim, hidden_dim, dropout) for i in range(num_layers - 2)]
        self.layern = AllConv(hidden_dim, input_edge_dim, output_node_dim, dropout)
        self.layers_store = torch.nn.Sequential(*self.layers)

    def forward(self, g, ts=None, groundtruth=None):
        with g.local_scope():
            nf = g.ndata['nf']
            x = self.layer0(g, nf)
            for layer in self.layers:
                x = layer(g, torch.cat([x, nf], dim=1)) + x
            x = self.layern(g, x)
            return None, None, x




class GAT(torch.nn.Module):
    def __init__(
        self,
        input_node_dim,
        input_edge_dim,
        hidden_dim,
        output_node_dim,
        num_layers,
        dropout,
    ):
        super().__init__()
        self.input_node_dim = input_node_dim
        self.input_edge_dim = input_edge_dim
        self.output_node_dim = output_node_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.dropout = dropout
        self.num_heads = 4
        
        self.layer_0 = dgl.nn.EGATConv(
            self.input_node_dim,
            self.input_edge_dim,
            self.hidden_dim,
            self.hidden_dim,
            self.num_heads,
        )
        
        self.layers = nn.ModuleList()
        for _ in range(self.num_layers-2):
            self.layers.append(
                dgl.nn.EGATConv(
                    self.hidden_dim,
                    self.hidden_dim,
                    self.hidden_dim,
                    self.hidden_dim,
                    self.num_heads,
                )
            )
        
        self.layer_n = dgl.nn.EGATConv(
            self.hidden_dim,
            self.hidden_dim,
            self.output_node_dim,
            self.hidden_dim,
            self.num_heads,
        )
        
        
    def forward(self, g, ts=None, groundtruth=None):
        with g.local_scope():
            nf = g.ndata['nf']
            ef = g.edata['ef']
            
            nf, ef = self.layer_0(g, nf, ef)
            nf = nf.mean(dim=1)
            ef = ef.mean(dim=1)
            
            for i in range(len(self.layers)):
                nf, ef = self.layers[i](g, nf, ef)
                nf = nf.mean(dim=1)
                ef = ef.mean(dim=1)
                
            nf, ef = self.layer_n(g, nf, ef)
            nf = nf.mean(dim=1)
            ef = ef.mean(dim=1)
        
        return None,None,nf



class GINE(torch.nn.Module):
    def __init__(
        self,
        input_node_dim,
        input_edge_dim,
        hidden_dim,
        output_node_dim,
        num_layers,
        dropout,
    ):
        super().__init__()
        self.input_node_dim = input_node_dim
        self.input_edge_dim = input_edge_dim
        self.output_node_dim = output_node_dim
        self.hidden_dim = hidden_dim
        self.num_layers = num_layers
        self.dropout = dropout

        # edge dim and node dim may not be the same, use MLP to make them match
        self.edge_to_node = MLP(self.input_edge_dim, self.hidden_dim, self.hidden_dim, self.hidden_dim, dropout=self.dropout)
        self.node_to_node = MLP(self.input_node_dim, self.hidden_dim, self.hidden_dim, self.hidden_dim, dropout=self.dropout)
    
        self.gine = torch.nn.ModuleList()
        self.gine.append(
            dgl.nn.GINEConv(MLP(self.hidden_dim, self.hidden_dim, self.hidden_dim, self.hidden_dim, dropout=self.dropout))
        )
        for _ in range(self.num_layers-1):
            self.gine.append(
                dgl.nn.GINEConv(MLP(self.hidden_dim, self.hidden_dim, self.hidden_dim, self.hidden_dim, dropout=self.dropout))
            )
        self.final_mlp = MLP(self.hidden_dim, self.hidden_dim, self.hidden_dim, self.output_node_dim, dropout=self.dropout)
    
    
    def forward(self, g, ts=None, groundtruth=None):
        with g.local_scope():
            node_feat = self.node_to_node(g.ndata['nf'])
            edge_feat = self.edge_to_node(g.edata['ef'])
            
            for i in range(len(self.gine)):
                node_feat = self.gine[i](g, node_feat, edge_feat)

            node_feat = self.final_mlp(node_feat)
        
        return None,None,node_feat