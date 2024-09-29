import torch
import dgl
import utils
import random
from . import process_data_hetero
from config import Config


def generate_homo_from_hetero_single(g_hetero):
    """dim of edge feature is 12"""
    na, nb = g_hetero.edges(etype='net_out', form='uv')
    ca, cb = g_hetero.edges(etype='cell_out', form='uv')

    # feature of net_out
    ne = torch.cat([
        torch.tensor([[0., 1., 0., 0., 0., 0., 0., 0., 0., 0.]]).expand(len(na), 10).to(g_hetero.device), 
        g_hetero.edges['net_out'].data['ef'][:,0:2]
    ], dim=1)

    # feature of cell_out
    ce = g_hetero.edges['cell_out'].data['ef'][:, 120:512].reshape(len(ca), 2*4, 49)
    ce = torch.cat([
        torch.tensor([[1., 0.]]).expand(len(ca), 2).to(g_hetero.device),
        torch.mean(ce, dim=2),
        torch.zeros(len(ca), 2).to(g_hetero.device)
    ], dim=1)

    g = dgl.graph((torch.cat([na, ca, nb, cb]), torch.cat([nb, cb, na, ca])))
    g.ndata['nf'] = g_hetero.ndata['nf']
    g.ndata['valid'] = g_hetero.ndata['valid']
    g.ndata['n_atslew'] = g_hetero.ndata['n_atslew']
    g.edata['ef'] = torch.cat([ne, ce, -ne, -ce])
    return g



def process_data_dict_homo(data_dict):
    data_dict_homo = {}
    for k, (g_hetero,_) in data_dict.items():
        g = generate_homo_from_hetero_single(g_hetero)
        if Config.model_type == "GAT":
            g = dgl.add_self_loop(g)
            print("model type = {}, add self loop".format(Config.model_type))
        data_dict_homo[k] = (g, {})

    return data_dict_homo



def process_data_homo(*args, **kwds):
    """
    for DeepGCNII baseline
    only predict AT and Slew
    will call data_hetero
    """
    data_train, data_test = process_data_hetero(*args, **kwds)
    data_train_homo = process_data_dict_homo(data_train)
    data_test_homo  = process_data_dict_homo(data_test)
    return data_train_homo, data_test_homo