from .block import Block
from .terminal import Terminal
from .fp_info import FPInfo
from typing import Tuple, Dict, Any, Union
import gymnasium as gym
import numpy as np
import torch
import math
import utils
from einops import rearrange, repeat
from copy import deepcopy
from matplotlib import pyplot as plt
from collections import OrderedDict
from tianshou.data import Batch
from collections import defaultdict
import os
import pandas as pd
plt.switch_backend('agg')


class RewardArgs:
    """
    Coefficients for reward calculation.
    """
    def __init__(self, 
                 reward_func:int,
                 reward_weight_hpwl:float, reward_weight_overlap:float, reward_weight_alignment:float, reward_weight_final_hpwl:float, 
                 ):
        self.reward_func = reward_func
        self.reward_weight_hpwl = reward_weight_hpwl
        self.reward_weight_overlap = reward_weight_overlap
        self.reward_weight_alignment = reward_weight_alignment
        self.reward_weight_final_hpwl = reward_weight_final_hpwl

    
    def __repr__(self) -> str:
        res = ""
        for key in self.__dict__:
            if not key.startswith("_"):
                res += "\t{:36} = {}\n".format(key, self.__dict__[key])
        return res


class PlaceEnv(gym.Env):
    def __init__(self, fp_info:FPInfo, 
                 overlap_ratio:float, along_boundary:bool, reward_args:RewardArgs, 
                 ratio_range:list[float,float], async_place:bool, 
                 device:torch.device, place_order_die_by_die:bool, input_next_block:int,
                 place_order_sorting_method:str,
                 graph:int, input_layer_sequence:bool, need_sequence_feature:bool,
                 need_alignment_mask:bool,
                 ):
        """
        Environment for floorplan.
        """
        super().__init__()
        # set variables
        self.fp_info = fp_info
        self.overlap_ratio = float(overlap_ratio)
        self.along_boundary = bool(along_boundary)
        self.ratio_range = ratio_range
        self.async_place = async_place
        self.device = device
        self.return_device = device
        self.place_order_die_by_die = place_order_die_by_die
        self.input_next_block = input_next_block
        self.reward_func = reward_args.reward_func
        self.place_order_sorting_method = place_order_sorting_method
        self.graph = graph
        self.input_layer_sequence = input_layer_sequence
        self.need_sequence_feature = need_sequence_feature
        self.empty_mask = torch.zeros((self.fp_info.x_grid_num, self.fp_info.y_grid_num)).to(device=self.return_device)
        self.need_alignment_mask = need_alignment_mask

        # define action space and observation space
        action_space = OrderedDict({"pos": gym.spaces.Discrete(fp_info.x_grid_num * fp_info.y_grid_num)})
        if ratio_range is not None:
            action_space["ratio"] = gym.spaces.Box(low=ratio_range[0], high=ratio_range[1], shape=(1,))
        if async_place:
            action_space["layer"] = gym.spaces.Discrete(self.fp_info.num_layer)

        self.action_space = gym.spaces.Dict(action_space)
        self.observation_space = None

        # set reward args
        self.reward_args = reward_args
        assert reward_args is not None, "[Error] reward_args is None."
        print(f"[INFO] reward_args: \n{self.reward_args}")

    
    def reset(self) -> Tuple[Dict, Dict]:
        # need to reset fp_info, blks and nets
        self.fp_info.reset()
        
        # placing order
        self.reset_place_order()

        self.reset_num_block_without_placing_order()

        if self.need_sequence_feature:
            self.reset_sequence_feature()
        
        if self.graph:
            self.reset_graph_data()

        self.last_metrics = defaultdict(lambda: 0)

        self.layer_curr_blk = 0 # the first block's die is 0

        # in next step, it is the current block
        next_block_moveable_idx = self.get_next_block_movable_idx(self.layer_curr_blk)
        assert next_block_moveable_idx is not None, "[Error] The first block is None on layer {}.".format(self.layer_curr_blk)
        next_block = self.fp_info.get_block_by_movable_idx(next_block_moveable_idx)
        self.update_num_block_without_placing_order(next_block)

        # one more block on this layer
        next_next_block_moveable_idx = self.get_next_next_block_movable_idx(self.layer_curr_blk)
        next_next_block = self.fp_info.get_block_by_movable_idx(next_next_block_moveable_idx) if next_next_block_moveable_idx is not None else None
        
        # reset layer sequence
        if self.input_layer_sequence:
            self.layer_sequence = np.ones(self.fp_info.movable_block_num, dtype=int) * -1 # -1 means this position is invalid
            self.layer_sequence_mask = np.ones(self.fp_info.movable_block_num, dtype=int) # 1 means this position is invalid
            self.layer_sequence_len = 0
            self.layer_sequence[0] = next_block.grid_z
            self.layer_sequence_mask[0] = 0
            self.layer_sequence_len += 1


        self.last_placed_block = {
            "area": -1,
            "height": -1,
            "placed": -1,
            "sequence": -1,
            "width": -1,
            "x": -1,
            "y": -1,
            "z": -1,
        }

        if self.need_alignment_mask:
            alignment_mask, binary_alignment_mask = self.get_alignment_mask(next_block, self.device)


        state = {
            "step": self.fp_info.placed_movable_block_num,
            "num_net": self.fp_info.net_num,
            "layer_idx": next_block.grid_z,
            "next_block_valid": 1,
            "ready_layers": self.get_ready_layers(),
            "num_blk_without_placing_order": self.num_block_without_placing_order.copy(),
            "last_placed_block": self.last_placed_block,

            "canvas": self.fp_info.canvas.to(device=self.return_device).clone(),
            "block": next_block,

            # mask
            "wiremask": self.get_wiremask(next_block, self.device).to(device=self.return_device),
            "position_mask": self.get_position_mask(next_block, self.along_boundary, device=self.device).to(device=self.return_device),
            "position_mask_loose": self.get_position_mask(next_block, False, self.overlap_ratio, device=self.device).to(device=self.return_device),
            "boundary_mask": self.get_boundary_mask(next_block, device=self.device).to(device=self.return_device),

            "wiremask_next": self.get_wiremask(next_next_block, self.device).to(device=self.return_device) if next_next_block is not None else self.empty_mask.clone(),
            "position_mask_next": self.get_position_mask(next_next_block, self.along_boundary, device=self.device).to(device=self.return_device) if next_next_block is not None else self.empty_mask.clone(),
            "grid_area_next": next_next_block.grid_area / self.get_mean_grid_area() if next_next_block is not None else 0,

        }

        # alignment mask
        if self.need_alignment_mask:
            state["alignment_mask"] = alignment_mask.to(device=self.return_device)
            state["binary_alignment_mask"] = binary_alignment_mask.to(device=self.return_device)
    

        # input_layer_sequence
        if self.input_layer_sequence:
            state["layer_sequence"] = self.layer_sequence
            state["layer_sequence_mask"] = self.layer_sequence_mask
            state["layer_sequence_len"] = self.layer_sequence_len
        
        # sequence feature
        if self.need_sequence_feature:
            state["sequence_feature"] = self.get_sequence_feature()
        
        # graph data
        if self.graph:
            state["graph_data"] = self.get_graph_data()
            state["graph_data"]["idx"] = torch.tensor(next_block.movable_idx)

        return state, {}


    
    def step(self, action: Union[OrderedDict,Batch]) -> tuple[Any, float, bool, bool, dict[str, Any]]:
        """state, reward, terminated, truncated, info"""
        # coordinates
        pos = action["pos"]
        x = pos // self.fp_info.y_grid_num
        y = pos % self.fp_info.y_grid_num

        # block ratio
        if self.ratio_range is not None:
            next_block_ratio = action["ratio"] if self.ratio_range is not None else None # next_block_ratio
            # next_block_ratio is in range [-1,1], use linear transformation to [low, high]
            next_block_ratio = (next_block_ratio + 1) / 2 * (self.ratio_range[1] - self.ratio_range[0]) + self.ratio_range[0]
            next_block_ratio = np.clip(next_block_ratio, self.ratio_range[0], self.ratio_range[1])
        

        # get block and pop, which is determined in last step
        curr_blk_mov_idx = self.get_next_block_movable_idx(self.layer_curr_blk, pop=True)
        assert curr_blk_mov_idx is not None, "[Error] Place a block, but the block is None (layer = {}).".format(self.layer_curr_blk)
        curr_blk = self.fp_info.get_block_by_movable_idx(curr_blk_mov_idx)

        # place block
        self.place_a_block(curr_blk, x, y)


        # done
        terminated = self.fp_info.is_all_placed()
        truncated = False


        # calculate reward
        action_z = action["layer"] if self.async_place else curr_blk.grid_z
        reward, info = self.calc_reward(terminated, x, y, action_z, self.fp_info.placed_movable_block_num, None)


        # layer_curr_blk, which is the layer to access next block
        if self.async_place:
            self.layer_curr_blk = action["layer"] # in asynchronous place, the layer of the next block is determined by the actor
            num_residual_blk = self.num_block_without_placing_order[self.layer_curr_blk]
            if num_residual_blk <= 0: # this die is invalid, change to another die
                layer_with_residual_blk = np.where(self.num_block_without_placing_order > 0)[0]
                if len(layer_with_residual_blk) > 0:
                    self.layer_curr_blk = layer_with_residual_blk[0]
                else:
                    self.layer_curr_blk = None
        else:
            self.layer_curr_blk = None # in synchronous place, the layer of the next block takes no effect

        # next block (in next step, this is the current block)
        next_block_moveable_idx = self.get_next_block_movable_idx(self.layer_curr_blk)
        next_block = self.fp_info.get_block_by_movable_idx(next_block_moveable_idx) if next_block_moveable_idx is not None else None
        if next_block is not None:
            self.layer_curr_blk = next_block.grid_z # next next block will be in the same die as the next block
            self.update_num_block_without_placing_order(next_block)

            # update layer sequence
            if self.input_layer_sequence:
                self.layer_sequence[self.fp_info.placed_movable_block_num] = next_block.grid_z
                self.layer_sequence_mask[self.fp_info.placed_movable_block_num] = 0
                self.layer_sequence_len += 1

            # update block ratio
            if self.ratio_range is not None:
                next_block.set_ratio(next_block_ratio, self.fp_info.x_grid_num, self.fp_info.y_grid_num)
        else:
            self.layer_curr_blk = None

        # next next block (one more block on this layer)
        next_next_block_layer = self.layer_curr_blk
        next_next_block_moveable_idx = self.get_next_next_block_movable_idx(next_next_block_layer)
        next_next_block = self.fp_info.get_block_by_movable_idx(next_next_block_moveable_idx) if next_next_block_moveable_idx is not None else None

        if next_block is not None:
            if self.need_alignment_mask:
                alignment_mask, binary_alignment_mask = self.get_alignment_mask(next_block, self.device)

        else:
            if self.need_alignment_mask:
                alignment_mask = self.empty_mask.clone()
                binary_alignment_mask = self.empty_mask.clone().to(dtype=torch.int32)

        obs_next = {
            "step": self.fp_info.placed_movable_block_num,
            "num_net": self.fp_info.net_num,
            "layer_idx": next_block.grid_z if next_block is not None else 0,
            "next_block_valid": 1 if next_block is not None else 0,
            "ready_layers": self.get_ready_layers(),
            "num_blk_without_placing_order": self.num_block_without_placing_order.copy(),
            "last_placed_block": self.last_placed_block,
            
            "canvas": self.fp_info.canvas.to(device=self.return_device).clone(),
            "block": next_block,

            # mask
            "wiremask": self.get_wiremask(next_block, self.device).to(device=self.return_device) if next_block is not None else self.empty_mask.clone(),
            "position_mask": self.get_position_mask(next_block, self.along_boundary, device=self.device).to(device=self.return_device) if next_block is not None else self.empty_mask.clone(),
            "position_mask_loose": self.get_position_mask(next_block, False, self.overlap_ratio, device=self.device).to(device=self.return_device) if next_block is not None else self.empty_mask.clone(),
            "boundary_mask": self.get_boundary_mask(next_block, self.device).to(device=self.return_device) if next_block is not None else self.empty_mask.clone(),

            "wiremask_next": self.get_wiremask(next_next_block, self.device).to(device=self.return_device) if next_next_block is not None else self.empty_mask.clone(),
            "position_mask_next": self.get_position_mask(next_next_block, self.along_boundary, device=self.device).to(device=self.return_device) if next_next_block is not None else self.empty_mask.clone(),
            "grid_area_next": next_next_block.grid_area / self.get_mean_grid_area() if next_next_block is not None else 0,
        }

        # alignment mask
        if self.need_alignment_mask:
            obs_next["alignment_mask"] = alignment_mask.to(device=self.return_device)
            obs_next["binary_alignment_mask"] = binary_alignment_mask.to(device=self.return_device)


        # input_layer_sequence
        if self.input_layer_sequence:
            obs_next["layer_sequence"] = self.layer_sequence
            obs_next["layer_sequence_mask"] = self.layer_sequence_mask
            obs_next["layer_sequence_len"] = self.layer_sequence_len
        
        # sequence feature
        if self.need_sequence_feature:
            obs_next["sequence_feature"] = self.get_sequence_feature()
        
        # graph data
        if self.graph:
            obs_next["graph_data"] = self.get_graph_data()
            obs_next["graph_data"]["idx"] = torch.tensor(next_block.movable_idx) if next_block is not None else torch.tensor(-1)

        return obs_next, reward, terminated, truncated, info
    

    def get_mean_grid_area(self) -> float:
        """
        Return the average block grid area.
        """
        if hasattr(self, "_mean_grid_area"):
            return self._mean_grid_area
        else:
            _mean_grid_area = 0.0
            for b in self.fp_info.block_info:
                _mean_grid_area += b.grid_area
            _mean_grid_area /= self.fp_info.block_num
            self._mean_grid_area = _mean_grid_area
            return _mean_grid_area


    def reset_graph_data(self):
        if not hasattr(self, "_has_init_graph_data"):
            self._has_init_graph_data = True
            n_preplaced = self.fp_info.preplaced_block_num
            n_movable = self.fp_info.movable_block_num
            # adjacency matrix
            s = n_preplaced
            e = n_preplaced + n_movable
            self.adj_mat_mov = self.fp_info.adjacency_matrix[s:e,s:e].clone().to(self.return_device)
            # set diagonal to 1
            self.adj_mat_mov[torch.eye(n_movable, dtype=bool)] = 1
            
            # node data
            # self.graph_x_init = torch.zeros(n_movable)
            # self.graph_y_init = torch.zeros(n_movable)
            # self.graph_z_init = torch.zeros(n_movable)
            # self.graph_w_init = torch.zeros(n_movable)
            # self.graph_h_init = torch.zeros(n_movable)
            # self.graph_area_init = torch.zeros(n_movable)
            # self.graph_placed_init = torch.zeros(n_movable)
            # if not self.async_place:
            #     self.graph_order_init = torch.zeros(n_movable)
            self.graph_x_init = np.zeros(n_movable)
            self.graph_y_init = np.zeros(n_movable)
            self.graph_z_init = np.zeros(n_movable)
            self.graph_w_init = np.zeros(n_movable)
            self.graph_h_init = np.zeros(n_movable)
            self.graph_area_init = np.zeros(n_movable)
            self.graph_placed_init = np.zeros(n_movable)
            if not self.async_place:
                self.graph_order_init = np.zeros(n_movable)

            for i in range(n_movable):
                block = self.fp_info.get_block_by_movable_idx(i)
                self.graph_x_init[i] = block.grid_x
                self.graph_y_init[i] = block.grid_y
                self.graph_z_init[i] = block.grid_z
                self.graph_w_init[i] = block.grid_w
                self.graph_h_init[i] = block.grid_h
                self.graph_area_init[i] = block.grid_area
                self.graph_placed_init[i] = block.placed
                if not self.async_place:
                    self.graph_order_init[i] = self.init_place_order.index(i)
        
        # reset
        # self.graph_x = self.graph_x_init.clone()
        # self.graph_y = self.graph_y_init.clone()
        # self.graph_z = self.graph_z_init.clone()
        # self.graph_w = self.graph_w_init.clone()
        # self.graph_h = self.graph_h_init.clone()
        # self.graph_area = self.graph_area_init.clone()
        # self.graph_placed = self.graph_placed_init.clone()
        # if not self.async_place:
        #     self.graph_order = self.graph_order_init.clone()
        self.graph_x = self.graph_x_init.copy()
        self.graph_y = self.graph_y_init.copy()
        self.graph_z = self.graph_z_init.copy()
        self.graph_w = self.graph_w_init.copy()
        self.graph_h = self.graph_h_init.copy()
        self.graph_area = self.graph_area_init.copy()
        self.graph_placed = self.graph_placed_init.copy()
        if not self.async_place:
            self.graph_order = self.graph_order_init.copy()

    def update_graph_data(self, block:Block):
        if block is None:
            return
        
        i = block.movable_idx
        self.graph_x[i] = block.grid_x
        self.graph_y[i] = block.grid_y
        self.graph_z[i] = block.grid_z
        self.graph_w[i] = block.grid_w
        self.graph_h[i] = block.grid_h
        self.graph_area[i] = block.grid_area
        self.graph_placed[i] = block.placed
    

    def get_graph_data(self) -> dict[str, torch.Tensor]:
        """
        Node feature of netlist graph.
        """
        graph_data = deepcopy(OrderedDict({
            "adj_mat_mov": self.adj_mat_mov,
            "x": self.graph_x,
            "y": self.graph_y,
            "z": self.graph_z,
            "w": self.graph_w,
            "h": self.graph_h,
            "area": self.graph_area,
            "placed": self.graph_placed,
        }))
        if not self.async_place:
            graph_data["order"] = self.graph_order
        return graph_data


    def set_hpwl_norm_coef(self, hpwl_norm_coef:float) -> None:
        print(f"[INFO] Set hpwl_norm_coef to {hpwl_norm_coef}.")
        self._hpwl_norm_coef = hpwl_norm_coef

    
    def reset_place_order(self):
        """
        async: place_order: List[List[int]].
        sync: place_order: List[int].
        The int in list is the movable index of block.
        """
        if hasattr(self, "init_place_order"):
            self.place_order = deepcopy(self.init_place_order)
            return 
        
        if self.place_order_sorting_method == "area":
            # synchronous place, place_order is a list of movable indices
            self.place_order = self.fp_info.get_unplaced_movable_block_movable_indices()
            self.place_order.sort(key=lambda movable_idx: self.fp_info.get_block_by_movable_idx(movable_idx).area, reverse=True)

            # asynchronous place, place_order is a list of list of movable indices
            if self.async_place:
                tmp = [ [] for _ in range(self.fp_info.num_layer) ]
                for idx in self.place_order:
                    block = self.fp_info.get_block_by_movable_idx(idx)
                    tmp[block.grid_z].append(idx)
                self.place_order = tmp
        

        elif self.place_order_sorting_method == "divide":
            blk_with_aln = [ blk_idx for blk_idx in self.fp_info.get_unplaced_movable_block_movable_indices() if self.fp_info.get_block_by_movable_idx(blk_idx).partner_indices.__len__() > 0 ]
            # sort with alignment area
            blk_with_aln.sort(key=lambda blk_idx: self.fp_info.get_block_by_movable_idx(blk_idx).alignment_area, reverse=True)

            blk_without_aln = [ blk_idx for blk_idx in self.fp_info.get_unplaced_movable_block_movable_indices() if self.fp_info.get_block_by_movable_idx(blk_idx).partner_indices.__len__() == 0 ]
            # sort with area
            blk_without_aln.sort(key=lambda blk_idx: self.fp_info.get_block_by_movable_idx(blk_idx).area, reverse=True)
            
            self.place_order = blk_with_aln + blk_without_aln

            if self.async_place:
                tmp = [ [] for _ in range(self.fp_info.num_layer) ]
                for idx in self.place_order:
                    block = self.fp_info.get_block_by_movable_idx(idx)
                    tmp[block.grid_z].append(idx)
                self.place_order = tmp

        else:            
            raise NotImplementedError(f"[Error] place_order_sorting_method {self.place_order_sorting_method} is not implemented.")


        # place_order_die_by_die
        if self.place_order_die_by_die:
            assert not self.async_place, "[Error] place_order_die_by_die is not supported in async_place."
            order = []
            for d in range(self.fp_info.num_layer):
                # get blocks in each die
                tmp = list(filter(lambda x: self.fp_info.get_block_by_movable_idx(x).grid_z == d, self.place_order))
                # sort by area
                tmp.sort(key=lambda x: self.fp_info.get_block_by_movable_idx(x).area, reverse=True)
                order.extend(tmp)
            self.place_order = order
        
        # move virtual block to the first
        self.virtual_blk = None
        for blk_idx in self.fp_info.get_unplaced_movable_block_movable_indices():
            blk = self.fp_info.get_block_by_movable_idx(blk_idx)
            if blk.virtual:
                self.virtual_blk = blk
                break
        if self.virtual_blk is not None:
            if self.async_place:
                virtual_blk_layer = self.virtual_blk.grid_z
                loc = self.place_order[virtual_blk_layer].index(self.virtual_blk.movable_idx)
                self.place_order[virtual_blk_layer].pop(loc)
                self.place_order[virtual_blk_layer].insert(0, self.virtual_blk.movable_idx)
            else:
                loc = self.place_order.index(self.virtual_blk.movable_idx)
                self.place_order.pop(loc)
                self.place_order.insert(0, self.virtual_blk.movable_idx)

            
        if self.async_place:
            # in each layer, unique
            for i in range(self.fp_info.num_layer):
                assert len(self.place_order[i]) == len(set(self.place_order[i]))
            # for arbitrary two layers, no intersection
            for i in range(self.fp_info.num_layer):
                for j in range(i+1, self.fp_info.num_layer):
                    assert len(set(self.place_order[i]) & set(self.place_order[j])) == 0
            # union of all layers is equal to all movable blocks
            assert set(sum(self.place_order, start=[])) == set(self.fp_info.get_unplaced_movable_block_movable_indices())
        else:
            # unique and equal
            assert len(self.place_order) == len(set(self.place_order))
            assert set(self.place_order) == set(self.fp_info.get_unplaced_movable_block_movable_indices())


        # copy to init_place_order
        self.init_place_order = deepcopy(self.place_order)
    

    def reset_sequence_feature(self):
        """sequence_feature: [num_die, C]"""
        if not hasattr(self, "init_seq_feat"):
            self.init_seq_feat = True

            # order of each die
            if self.async_place:
                max_len = max(map(len, self.place_order))
                sequence_each_die:list[list] = []
                for die in range(self.fp_info.num_layer):
                    seq = self.place_order[die] + [-1] * (max_len - len(self.place_order[die]))
                    sequence_each_die.append(seq)
            else:
                sequence_each_die = [ [] for _ in range(self.fp_info.num_layer) ]
                for idx in self.place_order:
                    block = self.fp_info.get_block_by_movable_idx(idx)
                    sequence_each_die[block.grid_z].append(idx)
                max_len = max(map(len, sequence_each_die))
                for die in range(self.fp_info.num_layer):
                    sequence_each_die[die] += [-1] * (max_len - len(sequence_each_die[die]))
            self.init_sequence_each_die = np.array(sequence_each_die) # [num_layer, max_len]

            # weight, height, area, placed
            self.init_width = np.full_like(self.init_sequence_each_die, -1)
            self.init_height = np.full_like(self.init_sequence_each_die, -1)
            self.init_area = np.full_like(self.init_sequence_each_die, -1)
            self.init_placed = np.full_like(self.init_sequence_each_die, -1)
            self.init_x = np.full_like(self.init_sequence_each_die, -1)
            self.init_y = np.full_like(self.init_sequence_each_die, -1)
            self.init_z = np.full_like(self.init_sequence_each_die, -1)
            for i in range(self.fp_info.num_layer):
                for j in range(len(self.init_sequence_each_die[i])):
                    if self.init_sequence_each_die[i][j] == -1: # not a valid block
                        continue
                    block = self.fp_info.get_block_by_movable_idx(self.init_sequence_each_die[i][j])
                    self.init_width[i,j] = block.grid_w
                    self.init_height[i,j] = block.grid_h
                    self.init_area[i,j] = block.grid_area
                    self.init_placed[i,j] = block.placed
                    self.init_x[i,j] = -1
                    self.init_y[i,j] = -1
                    self.init_z[i,j] = block.grid_z
        # reset
        self.sequence_each_die = self.init_sequence_each_die.copy()
        self.width = self.init_width.copy()
        self.height = self.init_height.copy()
        self.area = self.init_area.copy()
        self.placed = self.init_placed.copy()
        self.x = self.init_x.copy()
        self.y = self.init_y.copy()
        self.z = self.init_z.copy()
    

    def update_sequence_feature(self, block:Block):
        # update width, height, area, placed
        i,j = np.where(self.sequence_each_die == block.movable_idx)

        self.width[i,j] = block.grid_w
        self.height[i,j] = block.grid_h
        self.area[i,j] = block.grid_area
        self.placed[i,j] = block.placed
        self.x[i,j] = block.grid_x
        self.y[i,j] = block.grid_y
        self.z[i,j] = block.grid_z


    def get_sequence_feature(self) -> dict[str, np.ndarray]:
        """
        Return a dict with keys: sequence, width, height, area, placed, x, y, z. 
        Each feature is a numpy array with shape (num_layer, max_len).
        sequence: the place order of each die.
        """
        return deepcopy(OrderedDict({
            "sequence": self.sequence_each_die,
            "width": self.width,
            "height": self.height,
            "area": self.area,
            "placed": self.placed,
            "x": self.x,
            "y": self.y,
            "z": self.z,

            "last_placed_block": self.last_placed_block,
        }))
    
    def get_place_order_detailed_information(self) -> pd.DataFrame:
        """
        Detailed information for movable blocks.
        """
        df = pd.DataFrame()
        place_order = self.place_order if self.async_place else [self.place_order]
        for layer_idx in range(len(place_order)):
            for order, movable_idx in enumerate(place_order[layer_idx]):
                block = self.fp_info.get_block_by_movable_idx(movable_idx)
                line = pd.DataFrame([{
                    "idx": block.idx,
                    "movable_idx": movable_idx,
                    "name": block.name,
                    "order": order,
                    "w": block.w,
                    "h": block.h,
                    "z": block.z,
                    "area": block.area,
                    "alignment_blocks": [self.fp_info.get_module_by_full_idx(pid).name for pid in block.partner_indices],
                }])
                df = pd.concat([df, line], ignore_index=True)
        return df
    

    def get_ready_layers(self) -> np.ndarray:
        """Return a numpy array with shape (num_layer,). 1 if there are blocks whose placing order is not determined. Otherwise, 0."""
        return np.where(self.num_block_without_placing_order > 0, 1, 0)


    def reset_num_block_without_placing_order(self) -> None:
        """In each die, the placing order of each block is not determined."""
        if hasattr(self, "init_num_block_without_placing_order"):
            self.num_block_without_placing_order = self.init_num_block_without_placing_order.copy()
            return
        
        self.num_block_without_placing_order = np.zeros(self.fp_info.num_layer, dtype=int)
        for b in self.fp_info.block_info:
            if not b.preplaced:
                self.num_block_without_placing_order[b.grid_z] += 1
        self.init_num_block_without_placing_order = self.num_block_without_placing_order.copy()


    def update_num_block_without_placing_order(self, block:Block) -> None:
        """After selecting the next block to place, we will know the layer of the next block. So we need to update the num_block_without_placing_order."""
        layer = block.grid_z
        self.num_block_without_placing_order[layer] -= 1



    def get_next_block_movable_idx(self, layer:int, pop:bool=False) -> int:
        """
        If all blocks are placed, return None. Otherwise, return the next block movable index.
        layer: the layer of the next block, which only takes effect in async_place.
        pop: if True, pop the first block in the list, which means that the place order of this block is determined, even if it is not placed.
        """
        if self.async_place:
            # next layer is invalid, means all blocks are placed
            if layer is None:
                return None
            
            assert 0 <= layer < self.fp_info.num_layer, "[Error] layer out of range, got {}, but max_layer = {}".format(layer, self.fp_info.num_layer)

            # this layer has no block
            if len(self.place_order[layer]) == 0:
                return None
            
            if pop:
                return self.place_order[layer].pop(0)
            else:
                return self.place_order[layer][0]
        
        else:
            # next block does not exist
            if len(self.place_order) == 0:
                return None
            
            if pop:
                return self.place_order.pop(0)
            else:
                return self.place_order[0]



    def get_next_next_block_movable_idx(self, layer:int) -> int:
        """
        If available, return the next next block movable index. Otherwise, return None.
        layer: the layer of the next next block, which only takes effect in async_place.
        """
        if self.async_place:
            # next layer is invalid, means all blocks are placed
            if layer is None:
                return None
            
            assert 0 <= layer < self.fp_info.num_layer, "[Error] layer out of range, got {}, but max_layer = {}".format(layer, self.fp_info.num_layer)
            if len(self.place_order[layer]) <= 1: # must have at least 2 blocks
                return None
            return self.place_order[layer][1]
        
        else:
            if len(self.place_order) <= 1: # must have at least 2 blocks
                return None
            return self.place_order[1]



    def place_a_block(self, block:Block, x:int, y:int):
        assert not block.placed, "Block {} has been placed.".format(block.idx)

        # place block
        block.place(x, y)
        self.fp_info.placed_movable_block_num += 1

        # update net range
        for net in block.connected_nets:
            net.update(block)

        # update canvas
        self.fp_info.update_canvas(block)

        # update sequence feature
        if self.need_sequence_feature:
            self.update_sequence_feature(block)
        
        # update graph data
        if self.graph:
            self.update_graph_data(block)

        # update last placed block
        self.last_placed_block.update({
            "area": block.grid_area,
            "height": block.grid_h,
            "placed": block.placed,
            "sequence": block.movable_idx,
            "width": block.grid_w,
            "x": block.grid_x,
            "y": block.grid_y,
            "z": block.grid_z,
        })



    @property
    def curr_process(self) -> Tuple[int,int]:
        return self.fp_info.placed_movable_block_num, self.fp_info.movable_block_num
    

    def calc_reward(self, terminated:bool, action_x:int, action_y:int, action_z:int, step:int, tm:torch.Tensor) -> Tuple[float, dict]:
        # hpwl. Note that weight_hpwl is the weighted hpwl, which is involved in the reward calculation.
        hpwl, weight_hpwl, original_hpwl = self.fp_info.calc_hpwl()
        hpwl_norm_coef = self._hpwl_norm_coef if hasattr(self, "_hpwl_norm_coef") else 1
        hpwl_delta = (self.last_metrics["weight_hpwl"] - weight_hpwl) / hpwl_norm_coef

        # alignment
        if self.reward_args.reward_weight_alignment is not None:
            alignment_score = self.fp_info.calc_alignment_score()

        # overlap
        overlap = self.fp_info.get_overlap() # the lower the better

        # number of residual blocks in layer z
        num_residual_blk = self.num_block_without_placing_order[action_z]
        
                    
        if self.reward_func == 5: # Intermediate step, use difference of metrics. Final step, use final metrics. Recompute reward.
            if not terminated:
                last_weight_hpwl = self.last_metrics["weight_hpwl"]
                last_overlap = self.last_metrics["overlap"]

                reward = self.reward_args.reward_weight_hpwl * (last_weight_hpwl - weight_hpwl) / hpwl_norm_coef \
                    + self.reward_args.reward_weight_overlap * (last_overlap - overlap)
                                
                if self.reward_args.reward_weight_alignment is not None:
                    last_alignment_score = self.last_metrics["alignment"]
                    reward += self.reward_args.reward_weight_alignment * (alignment_score - last_alignment_score)
            

            else:
                reward = 0.0
                reward -= self.reward_args.reward_weight_overlap * overlap
                reward -= self.reward_args.reward_weight_final_hpwl * weight_hpwl / hpwl_norm_coef

                if self.reward_args.reward_weight_alignment is not None:
                    reward += self.reward_args.reward_weight_alignment * alignment_score
                

        else:
            raise NotImplementedError(f"[Error] reward_func: {self.reward_func} is not supported.")

        # calculate layer sum of first half sequence
        if step <= self.fp_info.movable_block_num // 2:
            layer_sum_first_half_seq = self.last_metrics["layer_sum_first_half_seq"] + action_z
        else:
            layer_sum_first_half_seq = self.last_metrics["layer_sum_first_half_seq"]

        # calculate whether next layer is valid
        next_layer_valid = 1 if num_residual_blk > 0 else 0
        next_layer_valid = self.last_metrics["next_layer_valid"] + next_layer_valid
        
        info = {
            "hpwl": hpwl,
            "hpwl_norm_coef": hpwl_norm_coef,
            "hpwl_delta": hpwl_delta,
            "weight_hpwl": weight_hpwl,
            "original_hpwl": original_hpwl,
            "overlap": overlap,
            "alignment": alignment_score if self.reward_args.reward_weight_alignment is not None else -1.0,

            "layer_sum_first_half_seq": layer_sum_first_half_seq,
            "next_layer_valid": next_layer_valid,
        }

        # update last_metrics
        self.last_metrics.update(info)

        return reward, info
    

    # @utils.record_time
    @torch.no_grad()
    def get_wiremask(self, tobe_placed_block:Block, device:torch.device=torch.device("cpu")) -> torch.Tensor:
        """wiremask.shape = (Nx, Ny)"""
        wiremask = torch.zeros((self.fp_info.x_grid_num, self.fp_info.y_grid_num), device=device)
        for net in tobe_placed_block.connected_nets:
            if net.num_placed_connector == 0:
                continue
            if 0 < net.x_min:
                x = torch.arange(0, net.x_min).to(device)
                x_rep = repeat(x, 'x -> x y', y=self.fp_info.y_grid_num)
                wiremask[x, :] += (net.x_min - x_rep) * net.get_net_weight()
            
            if net.x_max+1 < self.fp_info.x_grid_num:
                x = torch.arange(net.x_max+1, self.fp_info.x_grid_num).to(device)
                x_rep = repeat(x, 'x -> x y', y=self.fp_info.y_grid_num)
                wiremask[x, :] += (x_rep - net.x_max) * net.get_net_weight()

            if 0 < net.y_min:
                y = torch.arange(0, net.y_min).to(device)
                y_rep = repeat(y, 'y -> x y', x=self.fp_info.x_grid_num)
                wiremask[:, y] += (net.y_min - y_rep) * net.get_net_weight()

            if net.y_max+1 < self.fp_info.y_grid_num:
                y = torch.arange(net.y_max+1, self.fp_info.y_grid_num).to(device)
                y_rep = repeat(y, 'y -> x y', x=self.fp_info.x_grid_num)
                wiremask[:, y] += (y_rep - net.y_max) * net.get_net_weight()

        return wiremask



    # @utils.record_time
    @torch.no_grad()
    def get_position_mask(self, tobe_placed_block:Block, along_boundary:bool=True, overlap_ratio:float=0.0, device:torch.device=torch.device("cpu")) -> torch.Tensor:
        """
        A position mask which indicates available position for the tobe_placed_block.
        Tensor with shape (x_grid_num, y_grid_num).
        0: available.
        1: not available.
        """
        assert not tobe_placed_block.placed, "[Error] The block {} has been placed.".format(tobe_placed_block.idx)
        curr_layer = tobe_placed_block.grid_z
        w1 = tobe_placed_block.grid_w
        h1 = tobe_placed_block.grid_h


        overlap_w1 = round(w1 * overlap_ratio)
        overlap_h1 = round(h1 * overlap_ratio)

        if along_boundary is False:
            # all position are available
            position_mask = torch.zeros((self.fp_info.x_grid_num, self.fp_info.y_grid_num), device=device)
            
        else:
            # only the position along the boundary are available
            position_mask = torch.ones((self.fp_info.x_grid_num, self.fp_info.y_grid_num), device=device)

            for i in range(self.fp_info.block_num):
                if self.fp_info.block_info[i].virtual:
                    continue
                # a placed block, the same layer
                if self.fp_info.block_info[i].placed and self.fp_info.block_info[i].grid_z == curr_layer:
                    placed_block = self.fp_info.block_info[i]

                    x2 = placed_block.grid_x
                    y2 = placed_block.grid_y
                    w2 = placed_block.grid_w
                    h2 = placed_block.grid_h

                    overlap_w2 = round(w2 * overlap_ratio)
                    overlap_h2 = round(h2 * overlap_ratio)

                    min_overlap_w = min(overlap_w1, overlap_w2)
                    min_overlap_h = min(overlap_h1, overlap_h2)

                    # set the top region of placed block to available
                    # start_x = x2
                    start_x = x2 - w1 + 1
                    end_x   = x2 + w2
                    start_y = y2 + h2 - min_overlap_h
                    end_y   = y2 + h2

                    start_x = max(0, start_x)
                    start_y = max(0, start_y)
                    end_y  += 1
                    position_mask[start_x:end_x, start_y:end_y] = 0


                    # set the bottom region of placed block to available
                    # start_x = x2
                    start_x = x2 - w1 + 1
                    end_x   = x2 + w2
                    start_y = y2 - h1
                    end_y   = y2 - h1 + min_overlap_h

                    start_x = max(0, start_x)
                    start_y = max(0, start_y)
                    end_y  += 1
                    end_y = max(0, end_y)
                    position_mask[start_x:end_x, start_y:end_y] = 0


                    # set the left region of placed block to available
                    start_x = x2 - w1
                    end_x   = x2 - w1 + min_overlap_w
                    # start_y = y2
                    start_y = y2 - h1 + 1
                    end_y   = y2 + h2

                    start_x = max(0, start_x)
                    start_y = max(0, start_y)
                    end_x  += 1
                    end_x = max(0, end_x)
                    position_mask[start_x:end_x, start_y:end_y] = 0


                    # set the right region of placed block to available
                    start_x = x2 + w2 - min_overlap_w
                    end_x   = x2 + w2
                    # start_y = y2
                    start_y = y2 - h1 + 1
                    end_y   = y2 + h2

                    start_x = max(0, start_x)
                    start_y = max(0, start_y)
                    end_x  += 1
                    position_mask[start_x:end_x, start_y:end_y] = 0
            
            
            # set region to 1 due to boundary
            # left
            position_mask[0, :] = 0
            # right
            position_mask[self.fp_info.x_grid_num - w1, :] = 0
            # top
            position_mask[:, self.fp_info.y_grid_num - h1] = 0
            # bottom
            position_mask[:, 0] = 0

                    
        # set region to 1 due to placed blocks
        for i in range(self.fp_info.block_num):
            if self.fp_info.block_info[i].virtual:
                continue
            if self.fp_info.block_info[i].placed and self.fp_info.block_info[i].grid_z == curr_layer:
                placed_block = self.fp_info.block_info[i]

                x2 = placed_block.grid_x
                y2 = placed_block.grid_y
                w2 = placed_block.grid_w
                h2 = placed_block.grid_h

                overlap_w2 = round(w2 * overlap_ratio)
                overlap_h2 = round(h2 * overlap_ratio)

                min_overlap_w = min(overlap_w1, overlap_w2)
                min_overlap_h = min(overlap_h1, overlap_h2)

                # left
                start_x = x2 - w1 + min_overlap_w
                start_x += 1
                start_x = max(0, start_x)

                # bottom
                start_y = y2 - h1 + min_overlap_h
                start_y += 1
                start_y = max(0, start_y)

                # right
                end_x = x2 + w2 - min_overlap_w
                end_x = max(0, end_x)

                # top
                end_y = y2 + h2 - min_overlap_h
                end_y = max(0, end_y)

                position_mask[start_x:end_x, start_y:end_y] = 1


        # set region to 1 due to boundary
        start_x = self.fp_info.x_grid_num - w1
        start_x += 1
        start_y = self.fp_info.y_grid_num - h1
        start_y += 1
        position_mask[start_x:, :] = 1
        position_mask[:, start_y:] = 1

        return position_mask
    

    # @utils.record_time
    @torch.no_grad()
    def get_boundary_mask(self, block:Block, device:torch.device=torch.device("cpu")) -> torch.Tensor:
        """
        For the given block, ensure that it will not be placed outside the boundary. 
        Return a mask with shape (x_grid_num, y_grid_num). 
        1: not available. 0: available.
        """
        mask = torch.zeros((self.fp_info.x_grid_num, self.fp_info.y_grid_num), device=device)
        x_start = self.fp_info.x_grid_num - block.grid_w
        y_start = self.fp_info.y_grid_num - block.grid_h
        x_start += 1
        y_start += 1
        mask[x_start:, :] = 1
        mask[:, y_start:] = 1
        return mask
    
    @torch.no_grad()
    def get_alignment_mask(self, block_to_place:Block, device:torch.device=torch.device("cpu")) -> tuple[torch.Tensor, torch.IntTensor]:
        """
        Assume block B has some alignment partners, A1, A2, ..., An.
        We will generate an alignment mask and an binary alignment mask for each alignment pair (B, Ai).
        These masks will be merged to get the final alignment mask and binary alignment mask.
        Binary alignment mask is used to check whether the block B can be placed at the position (x, y).
        0 means the position is valid, 1 means the position is invalid.
        """
        alignment_masks = []
        binary_alignment_masks = []
        for pid in block_to_place.partner_indices:
            partner_block = self.fp_info.get_module_by_full_idx(pid)
            threshold = block_to_place.alignment_areas[partner_block.idx]
            alignment_mask = self._get_alignment_mask(block_to_place, partner_block, device)
            alignment_masks.append(alignment_mask)
            binary_alignment_mask = torch.where(alignment_mask >= threshold, 0, 1).to(dtype=torch.int32, device=device)
            binary_alignment_masks.append(binary_alignment_mask)
        
        if len(alignment_masks) == 0:
            single_alignment_mask = self._get_alignment_mask(block_to_place, None, device)
            single_binary_alignment_mask = torch.where(single_alignment_mask >= 0.0, 0, 1).to(dtype=torch.int32, device=device)
        else:
            single_alignment_mask = torch.stack(alignment_masks, dim=0).max(dim=0).values # union
            single_binary_alignment_mask = torch.stack(binary_alignment_masks, dim=0).prod(dim=0).to(dtype=torch.int32, device=device) # 0 * X = 0
        
        return single_alignment_mask, single_binary_alignment_mask

    
    @torch.no_grad()
    def _get_alignment_mask(self, block_to_place:Block, partner_block:Block, device:torch.device) -> torch.Tensor:
        """
        Return a mask with shape (x_grid_num, y_grid_num).
        Each element is the projection area of overlap between block_to_place and partner_block.
        """
        x1, y1, w1, h1 = block_to_place.grid_x, block_to_place.grid_y, block_to_place.grid_w, block_to_place.grid_h

        if partner_block is None or not partner_block.placed: # if not placed, we think every position is completely overlapped except the boundary
            mask = torch.zeros((self.fp_info.x_grid_num, self.fp_info.y_grid_num), device=device)
            x_start, y_start = 0, 0
            x_end, y_end = self.fp_info.x_grid_num - w1 + 1, self.fp_info.y_grid_num - h1 + 1
            mask[x_start:x_end, y_start:y_end] = block_to_place.grid_area
            return mask

        x2, y2, w2, h2 = partner_block.grid_x, partner_block.grid_y, partner_block.grid_w, partner_block.grid_h
        Nx, Ny = self.fp_info.x_grid_num, self.fp_info.y_grid_num

        x1_start = max(0, x2 - w1 + 1)
        x1_end = min(x2 + w2, Nx)
        y1_start = max(0, y2 - h1 + 1)
        y1_end = min(y2 + h2, Ny)

        overlap_fast = torch.zeros((Nx, Ny), device=device)
        x1 = torch.arange(x1_start, x1_end, device=device)
        y1 = torch.arange(y1_start, y1_end, device=device)
        x1, y1 = torch.meshgrid(x1, y1, indexing='ij')

        l1 = x1
        r1 = x1 + w1
        l2 = x2
        r2 = x2 + w2

        overlap_x = (torch.where(r1<r2, r1, r2) - torch.where(l1>l2, l1, l2)).clamp(min=0)

        l1 = y1
        r1 = y1 + h1
        l2 = y2
        r2 = y2 + h2

        overlap_y = (torch.where(r1<r2, r1, r2) - torch.where(l1>l2, l1, l2)).clamp(min=0)

        overlap_fast[x1,y1] = (overlap_x * overlap_y).to(overlap_fast.dtype)

        return overlap_fast

    
