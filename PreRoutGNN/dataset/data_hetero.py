import torch
import dgl
import random
import time
import utils
import os
from encoding import get_positional_encoding
from copy import deepcopy


def gen_topo(g_hetero):
    """
    generate topological levels from heterogeneous graph
    """
    na, nb = g_hetero.edges(etype='net_out', form='uv')
    ca, cb = g_hetero.edges(etype='cell_out', form='uv')
    g = dgl.graph((torch.cat([na, ca]).cpu(), torch.cat([nb, cb]).cpu()))
    topo = dgl.topological_nodes_generator(g)
    ret = [t.to(g_hetero.device) for t in topo]
    return ret


def generate_homo_from_hetero_single(g_hetero):
    """
    Generate homogeneous graph from heterogeneous graph.
    node_feature: keep the same with hetero
    net_out : 2 + 10, 10 is padding
    cell_out: 8 + 4 ,  4 is padding 
    """
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


def process_data_hetero(
        circuits_json_path:str, 
        predict_slew:bool, 
        predict_netdelay:bool,
        predict_celldelay:bool,
        data_split_json_path:str,
        num_level_freq_compents:int,
        sub_graph_size:int,
        num_pin_location_freq_compents:int,
        not_check_datasplit:bool,
        move_to_cuda_in_advance:bool,
        device:str,
        max_level:int,
        normalize_pin_location:bool,
        use_graph_autoencoder:bool,
        scale_capacitance:float,
    ):

    available_data = utils.load_json(circuits_json_path)
    if isinstance(data_split_json_path, str) and os.path.isfile(data_split_json_path):
        data_split = utils.load_json(data_split_json_path)
        if not not_check_datasplit:
            assert sorted(data_split["train"] + data_split["test"]) == sorted(available_data)
        train_data_keys = data_split["train"]
        test_data_keys = data_split["test"]
    else:
        train_data_keys = sorted(random.sample(available_data, int(0.8*len(available_data))))
        test_data_keys = sorted(list(set(available_data).difference(set(train_data_keys))))


    data = {}

    for circuit_name in available_data:
        g = dgl.load_graphs(circuit_name)[0][0]
        g.ndata['valid'] = torch.ones(g.num_nodes())
        
        if not predict_slew: del g.ndata['n_slews']
        if not predict_netdelay: del g.ndata['n_net_delays']
        if not predict_celldelay: del g.edges['cell_out'].data['e_cell_delays']
        
        # scale capacitance
        other_part = g.ndata['nf'][:, 0:6]
        capacitance = g.ndata['nf'][:, 6:10] * scale_capacitance
        g.ndata['nf'] = torch.cat([other_part, capacitance], dim=1)
        
        
        # normalize pin location through min_max normalization
        if normalize_pin_location:
            part1 = g.ndata['nf'][:, 0:2]
            loc   = g.ndata['nf'][:, 2:6]
            part2 = g.ndata['nf'][:, 6:10]
            loc = (loc - loc.min(dim=0)[0]) / (loc.max(dim=0)[0] - loc.min(dim=0)[0] + 1e-8)
            g.ndata['nf'] = torch.cat([part1, loc, part2], dim=1)


        if move_to_cuda_in_advance:
            g = g.to(device)


        if predict_netdelay: 
            g.ndata['n_net_delays_log'] = torch.log(0.0001 + g.ndata['n_net_delays']) + 7.6
        
        invalid_nodes = torch.abs(g.ndata['n_ats']) > 1e20
        g.ndata['n_ats'][invalid_nodes] = 0
        if predict_slew: 
            g.ndata['n_slews'] = torch.log(0.0001 + g.ndata['n_slews']) + 3

        

        if predict_slew:
            g.ndata['n_atslew'] = torch.cat([
                g.ndata['n_ats'],
                g.ndata['n_slews'],
            ], dim=1)
        else:
            g.ndata['n_atslew'] = g.ndata['n_ats']


        # dtype is float32
        g.edges['cell_out'].data['ef'] = g.edges['cell_out'].data['ef'].type(torch.float32)
        if predict_celldelay: g.edges['cell_out'].data['e_cell_delays'] = g.edges['cell_out'].data['e_cell_delays'].type(torch.float32)

        topo = gen_topo(g)
        
        if num_level_freq_compents > 0:
            levels = torch.zeros(g.num_nodes()).to(g.ndata['nf'].device)
            for level, nodes in enumerate(topo):
                levels[nodes.long()] = level
            level_encoding = get_positional_encoding(levels, num_level_freq_compents, max_level)
            level_encoding = torch.cat([levels.unsqueeze(-1), level_encoding], dim=1)
            g.ndata['nf'] = torch.cat([g.ndata['nf'], level_encoding], dim=1)
        
        
        if num_pin_location_freq_compents > 0:
            # pin
            pin_location_encoding = torch.cat([
                get_positional_encoding(g.ndata['nf'][:,i], num_pin_location_freq_compents) for i in range(2,6)
            ], dim=1)
            g.ndata['nf'] = torch.cat([g.ndata['nf'], pin_location_encoding], dim=1)
            
            # net_out
            pin_location_encoding = torch.cat([
                get_positional_encoding(g.edges['net_out'].data['ef'][:,i], num_pin_location_freq_compents) for i in range(0,2)
            ], dim=1)
            g.edges['net_out'].data['ef'] = torch.cat([g.edges['net_out'].data['ef'], pin_location_encoding], dim=1)
            
            # net in
            pin_location_encoding = torch.cat([
                get_positional_encoding(g.edges['net_in'].data['ef'][:,i], num_pin_location_freq_compents) for i in range(0,2)
            ], dim=1)
            g.edges['net_in'].data['ef'] = torch.cat([g.edges['net_in'].data['ef'], pin_location_encoding], dim=1)


        ts = {
            'input_nodes': (g.ndata['nf'][:, 1] < 0.5).nonzero().flatten().type(torch.int32),
            'output_nodes': (g.ndata['nf'][:, 1] > 0.5).nonzero().flatten().type(torch.int32),
            'output_nodes_nonpi': torch.logical_and(g.ndata['nf'][:, 1] > 0.5, g.ndata['nf'][:, 0] < 0.5).nonzero().flatten().type(torch.int32),
            'pi_nodes': torch.logical_and(g.ndata['nf'][:, 1] > 0.5, g.ndata['nf'][:, 0] > 0.5).nonzero().flatten().type(torch.int32),
            'po_nodes': torch.logical_and(g.ndata['nf'][:, 1] < 0.5, g.ndata['nf'][:, 0] > 0.5).nonzero().flatten().type(torch.int32),
            'endpoints': (g.ndata['n_is_timing_endpt'] > 0.5).nonzero().flatten().type(torch.int32),
            'topo': topo,
            'name': circuit_name,
            'valid': g.nodes().long(), # no partition, all nodes are valid
            'homo': generate_homo_from_hetero_single(g) if use_graph_autoencoder else None,
        }
                
        data[circuit_name] = g, ts


    data_train = {k: t for k, t in data.items() if k in train_data_keys}
    data_test = {k: t for k, t in data.items() if k in test_data_keys}

    
    if sub_graph_size > 0:
        data_train_partition = {}
        
        # only do partition on training set
        for circuit_name, (g, ts) in data_train.items():
            original_topo = ts["topo"]
            num_levels = len(original_topo)
            
            if g.num_nodes() <= sub_graph_size:
                data_train_partition[circuit_name] = g, ts
                
            else:
                i = 0
                start_level_id = 0
                while start_level_id < num_levels:
                    # fix start_level_id, increase block_level_size
                    block_level_size = 2
                    valid_nodes = torch.cat(original_topo[start_level_id: start_level_id+block_level_size], dim=0)
                    while valid_nodes.shape[0] < sub_graph_size and start_level_id+block_level_size+2 <= num_levels:
                        block_level_size += 2
                        valid_nodes = torch.cat(original_topo[start_level_id: start_level_id+block_level_size], dim=0)
                        
                    g_copy = deepcopy(g)
                    g_copy.ndata['valid'] = torch.zeros(g_copy.num_nodes()).to(g.device)
                    g_copy.ndata['valid'][valid_nodes.long()] = 1.0
                    
                    # add padding nodes
                    sub_graph_nodes = torch.cat(original_topo[max(start_level_id-4, 0) : start_level_id+block_level_size+4], dim=0)
                    sub_graph = dgl.node_subgraph(g_copy, sub_graph_nodes)
                    
                    small_topo = gen_topo(sub_graph)
                    sub_graph_name = f"{circuit_name}-{i}"
                    sub_ts = {
                        'input_nodes': (sub_graph.ndata['nf'][:, 1] < 0.5).nonzero().flatten().type(torch.int32),
                        'output_nodes': (sub_graph.ndata['nf'][:, 1] > 0.5).nonzero().flatten().type(torch.int32),
                        'output_nodes_nonpi': torch.logical_and(sub_graph.ndata['nf'][:, 1] > 0.5, sub_graph.ndata['nf'][:, 0] < 0.5).nonzero().flatten().type(torch.int32),
                        'pi_nodes': torch.logical_and(sub_graph.ndata['nf'][:, 1] > 0.5, sub_graph.ndata['nf'][:, 0] > 0.5).nonzero().flatten().type(torch.int32),
                        'po_nodes': torch.logical_and(sub_graph.ndata['nf'][:, 1] < 0.5, sub_graph.ndata['nf'][:, 0] > 0.5).nonzero().flatten().type(torch.int32),
                        'endpoints': (sub_graph.ndata['n_is_timing_endpt'] > 0.5).nonzero().flatten().type(torch.int32),
                        'topo': small_topo,
                        'name': sub_graph_name,
                        'valid': torch.where(sub_graph.ndata['valid'] > 0.5)[0].to(sub_graph.device).long(),
                        'homo': generate_homo_from_hetero_single(sub_graph) if use_graph_autoencoder else None,
                    }
                    
                    data_train_partition[sub_graph_name] = sub_graph, sub_ts
                    i += 1
                    start_level_id += block_level_size
    
    
    if sub_graph_size > 0: 
        data_train = data_train_partition
    
    
    # for k in sorted(data_train.keys()):
    #     print("Training {:50} num_nodes = {:6d}, num_levels = {}".format(k, data_train[k][0].num_nodes(), len(data_train[k][1]['topo'])))
    # for k in sorted(data_test.keys()):
    #     print("Test     {:50} num_nodes = {:6d}, num_levels = {}".format(k, data_test[k][0].num_nodes(),  len(data_test[k][1]['topo'])))
    
    print(f"Number of graphs in training: {len(data_train)}")
    print(f"Number of graphs in testing : {len(data_test)}")
            

    return data_train, data_test