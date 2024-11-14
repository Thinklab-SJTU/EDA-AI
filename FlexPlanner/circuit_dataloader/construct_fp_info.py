from .parser import parse_blk_tml, map_tml, parse_net
from .construct_partner import construct_partner_blk
from .construct_layer import assign_layer
from .construct_pre_placed_module import construct_preplaced_modules
import torch
import fp_env
import pandas as pd
from typing import Tuple
from copy import deepcopy
from collections import defaultdict


def construct_fp_info_func(circuit:str, area_util:float, num_grid_x:int, num_grid_y:int, num_alignment:int, 
                           alignment_rate:float, alignment_sort:str, num_preplaced_module:int, add_virtual_block:bool, num_layer:int) -> Tuple[fp_env.FPInfo, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    alignment_rate = 1.0 if alignment_rate is None else alignment_rate
    # read block and terminal (w,h,virtual)
    blk_wh_dict, tml_xy_dict, outline_width, outline_height = parse_blk_tml(circuit, area_util)
    tml_xy_dict = map_tml(tml_xy_dict, outline_width, outline_height)

    # assign layer (w,h,virtual,z)
    first_blk_wh = next(iter(blk_wh_dict.values()))
    if "z" not in first_blk_wh.keys():
        print("[INFO] Constructing layer information")
        blk_wh_dict = assign_layer(blk_wh_dict, num_layer)

    # construct preplaced modules (w,h,virtual,z,preplaced,x,y)
    first_blk_wh = next(iter(blk_wh_dict.values()))
    if "preplaced" not in first_blk_wh.keys():
        print("[INFO] Constructing preplaced modules")
        blk_wh_dict = construct_preplaced_modules(num_preplaced_module, blk_wh_dict, outline_height)

    # add virtual block
    if add_virtual_block:
        print("[INFO] Adding virtual block")
        blk_wh_dict["virtual_block"] = {
            "w": 1,
            "h": 1,
            "virtual": True,
            "z": 0,
            "preplaced": False,
            "x": 0,
            "y": 0,
        }


    # all blocks, preplaced_blocks + movable_blocks
    preplaced_blocks = []
    movable_blocks = []
    for blk_name, blk_info in blk_wh_dict.items():
        if blk_info['preplaced']: # PPM
            preplaced_blocks.append(fp_env.Block(blk_info['x'], blk_info['y'], blk_info['z'], blk_info['w'], blk_info['h'], blk_name, True, blk_info['virtual']))
        else: # movable or virtual
            movable_blocks.append(fp_env.Block(0, 0, blk_info['z'], blk_info['w'], blk_info['h'], blk_name, False, blk_info['virtual']))

    block_info = preplaced_blocks + movable_blocks


    # terminal
    terminal_info = []
    for terminal_name in tml_xy_dict.keys():
        terminal_xy = tml_xy_dict[terminal_name]
        terminal_info.append(fp_env.Terminal(terminal_xy['x'], terminal_xy['y'], 0, terminal_name))


    # discretize block and terminal
    block_info, terminal_info, grid_width, grid_height = fp_env.discretize(block_info, terminal_info, num_grid_x, num_grid_y, outline_width, outline_height)


    # read nets and construct net_info, nets is List[List[str]]
    # also construct adjacency matrix
    nets_str = parse_net(circuit)
    name2obj = {obj.name: obj for obj in block_info + terminal_info}
    net_info = []
    net_weight = 1.0
    for net_connectors_str in nets_str:
        connector_list = [name2obj[connector_name] for connector_name in net_connectors_str]
        net = fp_env.Net(connector_list, net_weight)
        net_info.append(net)


    # load to fp_info
    fp_info = fp_env.FPInfo(block_info, terminal_info, net_info, outline_width, outline_height, num_grid_x, num_grid_y)
    episode_len = fp_info.movable_block_num

    # construct adjacency matrix
    adjacency_matrix = torch.zeros(fp_info.block_num + fp_info.termimal_num, fp_info.block_num + fp_info.termimal_num)
    for net_idx, net_connectors_str in enumerate(nets_str):
        for i in range(len(net_connectors_str)):
            for j in range(i+1, len(net_connectors_str)):
                adjacency_matrix[name2obj[net_connectors_str[i]].idx, name2obj[net_connectors_str[j]].idx] = fp_info.net_info[net_idx].get_net_weight()
                adjacency_matrix[name2obj[net_connectors_str[j]].idx, name2obj[net_connectors_str[i]].idx] = fp_info.net_info[net_idx].get_net_weight()
            # add connector to net
            fp_info.net_info[net_idx].add_connector(name2obj[net_connectors_str[i]])
    fp_info.set_adjacency_matrix(adjacency_matrix)

    # construct alignment partner, and set to fp_info
    fp_info.set_alignment_sort(alignment_sort)
    df_partner = construct_partner_blk(fp_info, num_alignment, alignment_sort)
    for row in df_partner.itertuples():
        blk0_name, blk1_name = row.blk0, row.blk1
        blk0, blk1 = name2obj[blk0_name], name2obj[blk1_name]
        blk0_idx, blk1_idx = blk0.idx, blk1.idx
        alignment_area = min(blk0.area, blk1.area) * alignment_rate
        fp_info.set_partner(blk0_idx, blk1_idx, alignment_area)
    
    for blk_name, aln_group in fp_info.name2alignment_group.items():
        df_partner.loc[df_partner['blk0'] == blk_name, 'alignment_group'] = aln_group
    df_partner.sort_values(by=['alignment_group', 'blk0', 'blk1'], inplace=True)
    df_partner["alignment_group"] = df_partner["alignment_group"].astype(int)

    return fp_info, df_partner