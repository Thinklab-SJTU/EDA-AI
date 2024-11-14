from .block import Block
from .net import Net
from .terminal import Terminal
from typing import List, Tuple, Union
import torch
import pandas as pd
from collections import Counter, defaultdict
import numpy as np


def discretize(block_info:List[Block], terminal_info:List[Terminal], x_grid_num:int, y_grid_num:int, original_outline_width:float, original_outline_height:float):
        """discretize the block x,y,w,h to grid position."""
        print('[INFO] discretize block to {} x {}'.format(x_grid_num, y_grid_num))

        grid_width  = original_outline_width  / x_grid_num
        grid_height = original_outline_height / y_grid_num

        canvas = defaultdict(lambda: torch.zeros(x_grid_num, y_grid_num).to(torch.int))

        for block in block_info:
            block.set_grid_wh(grid_width, grid_height)
            block.set_grid_xy(grid_width, grid_height, x_grid_num, y_grid_num)
            if block.preplaced:
                while True:
                    canvas[block.grid_z][block.grid_x:block.grid_x+block.grid_w, block.grid_y:block.grid_y+block.grid_h] += 1
                    if (canvas[block.grid_z] > 1).sum() > 0:
                        canvas[block.grid_z][block.grid_x:block.grid_x+block.grid_w, block.grid_y:block.grid_y+block.grid_h] -= 1
                        block.grid_x += 1
                    else:
                        break
                    
        
        print('[INFO] discretize terminal to {} x {}'.format(x_grid_num, y_grid_num))
        for terminal in terminal_info:
            terminal.set_grid_xy(grid_width, grid_height, x_grid_num, y_grid_num)
        

        return block_info, terminal_info, grid_width, grid_height


class FPInfo:
    def __init__(self, block_info:List[Block], terminal_info:List[Terminal], net_info:List[Net], \
                 original_outline_width:float, original_outline_height:float, x_grid_num:int, y_grid_num:int):
        """original_outline_width and original_outline_height are the original width and height of the chip. They are used to discretize the block position to grid position."""
        self.block_info = block_info
        self.terminal_info = terminal_info
        self.net_info = net_info

        self.original_outline_width = original_outline_width
        self.original_outline_height = original_outline_height

        self.x_grid_num = x_grid_num
        self.y_grid_num = y_grid_num

        self.grid_width = original_outline_width / x_grid_num
        self.grid_height = original_outline_height / y_grid_num
        print('[INFO] grid_width: {:.3f}'.format(self.grid_width))
        print('[INFO] grid_height: {:.3f}'.format(self.grid_height))
        print('[INFO] grid_area: {:.3f}'.format(self.grid_width * self.grid_height))

        self.block_num = len(block_info)
        self.termimal_num = len(terminal_info)
        self.net_num = len(net_info)


        assert self.block_num > 0, "[ERROR] block_num should be larger than 0, but got {}".format(self.block_num)

        # for all partner pairs
        self.partner_idx2indices: dict[int, list[int]] = defaultdict(list)

        # list of add modules
        self.all_modules: List[Union[Block, Terminal]] = block_info + terminal_info

        # set index for block
        self.preplaced_block_indices = [] # full index for preplaced block
        self.movable_block_indices = [] # full index for movable block
        preplaced_block_idx = 0
        movable_block_idx = 0

        for full_idx, block in enumerate(self.block_info):
            if block.preplaced:
                self.preplaced_block_indices.append(full_idx)
                block.set_idx(full_idx, preplaced_block_idx)
                preplaced_block_idx += 1
            else:
                self.movable_block_indices.append(full_idx)
                block.set_idx(full_idx, movable_block_idx)
                movable_block_idx += 1
        
        self.preplaced_block_num = len(self.preplaced_block_indices)
        self.movable_block_num = len(self.movable_block_indices)
        print('[INFO] preplaced_block_num: {}'.format(self.preplaced_block_num))
        print('[INFO] movable_block_num: {}'.format(self.movable_block_num))

        # set terminal full idx and terminal_idx
        for tml_idx, tml in enumerate(self.terminal_info):
            full_idx = tml_idx + self.block_num
            tml.set_idx(full_idx, tml_idx)
        print('[INFO] terminal_num: {}'.format(self.termimal_num))


        # num_layer
        self.num_layer = 0
        for connector in self.block_info + self.terminal_info:
            self.num_layer = max(self.num_layer, connector.z + 1)
        print('[INFO] num_layer: {}'.format(self.num_layer))
        print('[INFO] net_num: {}'.format(self.net_num))

    def set_alignment_sort(self, alignment_sort:str):
        self.alignment_sort = alignment_sort
    
    def set_adjacency_matrix(self, adjacency_matrix:torch.Tensor):
        """set adjacency matrix for the net_info."""
        self.adjacency_matrix = adjacency_matrix
        print('[INFO] adjacency_matrix: {}'.format(self.adjacency_matrix.shape))
        print('[INFO] mean degree: {}'.format(self.adjacency_matrix.sum() / self.adjacency_matrix.shape[0]))
    

    @torch.no_grad()
    def get_overlap(self, norm:bool=True) -> int:
        overlap = (self.canvas - 1).clamp_min_(0).sum().item()
        if norm:
            overlap /= (self.canvas.shape[1] * self.canvas.shape[2])
        return overlap


    def reset_canvas(self):
        """
        canvas.shape is (num_layer, x_grid_num, y_grid_num)
        """
        self.canvas = torch.zeros(self.num_layer, self.x_grid_num, self.y_grid_num)
        for block in self.block_info:
            if block.preplaced and block.placed:
                self.canvas[block.grid_z, block.grid_x:block.grid_x+block.grid_w, block.grid_y:block.grid_y+block.grid_h] += 1
    

    def update_canvas(self, block:Block):
        """update canvas with the new placed block_idx."""
        if block.virtual:
            return
        self.canvas[block.grid_z, block.grid_x:block.grid_x+block.grid_w, block.grid_y:block.grid_y+block.grid_h] += 1


    def reset(self):
        """
        set placed to False for movable blocks.
        set net range.
        """
        for block in self.block_info:
            block.reset()
        
        for net in self.net_info:
            net.reset()
        
        self.placed_movable_block_num = 0

        self.reset_canvas()



    def calc_hpwl(self) -> Tuple[float, float, float]:
        """calculate grid_hpwl, weighted_grid_hpwl, original_hpwl."""
        grid_hpwl = 0
        weighted_grid_hpwl = 0
        original_hpwl = 0

        for net in self.net_info:
            tmp = net.calc_hpwl()
            grid_hpwl += tmp
            weighted_grid_hpwl += tmp * net.get_net_weight()
            x_stride, y_stride = net.calc_stride()
            original_hpwl += x_stride * self.grid_width + y_stride * self.grid_height

        return grid_hpwl, weighted_grid_hpwl, original_hpwl
    

    def calc_original_hpwl(self) -> float:
        """calculate HPWL for all nets."""
        hpwl = 0
        for net in self.net_info:
            x_stride, y_stride = net.calc_stride()
            hpwl += x_stride * self.grid_width + y_stride * self.grid_height
        return hpwl


    def get_block_by_movable_idx(self, movable_idx:int) -> Block:
        assert 0 <= movable_idx < self.movable_block_num, "movable_idx: {}, movable_block_num: {}".format(movable_idx, self.movable_block_num)
        return self.block_info[self.movable_block_indices[movable_idx]]
    

    def get_block_by_preplaced_idx(self, preplaced_idx:int) -> Block:
        assert 0 <= preplaced_idx < self.preplaced_block_num, "preplaced_idx: {}, preplaced_block_num: {}".format(preplaced_idx, self.preplaced_block_num)
        return self.block_info[self.preplaced_block_indices[preplaced_idx]]
    

    def get_terminal_by_terminal_idx(self, terminal_idx:int) -> Terminal:
        """use the terminal_idx, instead of full idx"""
        assert 0 <= terminal_idx < self.termimal_num, "terminal_idx: {}, termimal_num: {}".format(terminal_idx, self.termimal_num)
        return self.terminal_info[terminal_idx]
    

    def get_unplaced_movable_block_movable_indices(self) -> List[int]:
        """return the movable indices, which are unplaced."""
        return [block.movable_idx for block in self.block_info if not block.placed and not block.preplaced]
    

    def get_module_by_full_idx(self, full_idx:int) -> Union[Block, Terminal]:
        """
        Get module by full idx.
        Module is either block or terminal.
        """
        assert 0 <= full_idx < self.block_num + self.termimal_num, "full_idx: {}, block_num: {}, termimal_num: {}".format(full_idx, self.block_num, self.termimal_num)
        return self.all_modules[full_idx]
    

    def is_all_placed(self) -> bool:
        # return all([block.placed for block in self.block_info])
        return self.placed_movable_block_num == self.movable_block_num
    

    def set_partner(self, curr_idx:int, part_idx:int, alignment_area:float):
        """
        set partner for curr_idx and part_idx.
        idx is the full idx for block.
        alignment_area is the original area for alignment.
        """
        assert 0 <= curr_idx < self.block_num, "curr_idx={} is out of range [0, {})".format(curr_idx, self.block_num)
        assert 0 <= part_idx < self.block_num, "part_idx={} is out of range [0, {})".format(part_idx, self.block_num)
        assert curr_idx != part_idx
        assert self.block_info[curr_idx].grid_z != self.block_info[part_idx].grid_z, "curr_idx: {}, layer = {}, part_idx: {}, layer = {}".format(curr_idx, self.block_info[curr_idx].grid_z, part_idx, self.block_info[part_idx].grid_z)
        self.partner_idx2indices[curr_idx].append(part_idx)
        self.partner_idx2indices[part_idx].append(curr_idx)

        # discretize the alignment area
        original_alignment_area = alignment_area
        alignment_area = alignment_area // (self.grid_width * self.grid_height)
        alignment_area = max(1, alignment_area)

        self.block_info[curr_idx].set_partner(part_idx, alignment_area, original_alignment_area)
        self.block_info[part_idx].set_partner(curr_idx, alignment_area, original_alignment_area)

        # print('[INFO] set partner: {} and {}, alignment_area: {}'.format(curr_idx, part_idx, alignment_area))

    
    def calc_alignment_score(self) -> float:
        """
        calculate the alignment area for all partners.
        return the average alignment rate.
        """
        if len(self.partner_idx2indices) == 0:
            return 0.0
        
        def calc_overlap_1d(x1:int, x2:int, w1:int, w2:int) -> int:
            """calculate the overlap area for 1D."""
            left, right = max(x1, x2), min(x1 + w1, x2 + w2)
            return max(0, right - left)
        
        sum_alignment_score = 0.0
        alignment_num = 0

        for blk_idx in self.partner_idx2indices.keys():
            block = self.get_module_by_full_idx(blk_idx)
            if not block.placed or len(block.partner_indices) == 0:
                continue
            
            # for this block, one aligns to multiple partners
            alignment_area = 0
            required_alignment_area = 0
            for partner_idx in block.partner_indices:
                partner_block = self.get_module_by_full_idx(partner_idx)
                if partner_block.placed:
                    alignment_area += calc_overlap_1d(block.grid_x, partner_block.grid_x, block.grid_w, partner_block.grid_w) * \
                                    calc_overlap_1d(block.grid_y, partner_block.grid_y, block.grid_h, partner_block.grid_h)
                    required_alignment_area += partner_block.grid_area
            
            required_alignment_area = min(required_alignment_area, block.grid_area)
            if required_alignment_area > 0: # at least one partner is placed
                alignment_score = np.clip(alignment_area / required_alignment_area, 0.0, 1.0)
                sum_alignment_score += alignment_score
                alignment_num += 1

        
        return sum_alignment_score / alignment_num if alignment_num > 0 else 0.0
    

    @property
    def name2alignment_group(self) -> dict[str, int]:
        if hasattr(self, '_name2alignment_group'):
            return self._name2alignment_group

        blks = sorted(list(set( blk.name for blk in self.block_info if blk.partner_indices.__len__() > 0 )))
        blk2int = {blk: i for i, blk in enumerate(blks)}
        int2blk = {i: blk for i, blk in enumerate(blks)}

        adj = torch.zeros(len(blks), len(blks))
        for blk in self.block_info:
            if blk.partner_indices.__len__() > 0:
                for partner_idx in blk.partner_indices:
                    partner_blk = self.get_module_by_full_idx(partner_idx)
                    partner_name = partner_blk.name

                    adj[blk2int[blk.name], blk2int[partner_name]] = 1
                    adj[blk2int[partner_name], blk2int[blk.name]] = 1
        
        # BFS this graph
        visited = set()
        name2alignment_group = {}
        group_idx = 0
        for blk in blks:
            if blk in visited:
                continue
            queue = [blk]
            while queue:
                curr_blk = queue.pop(0)
                if curr_blk in visited:
                    continue
                visited.add(curr_blk)
                name2alignment_group[curr_blk] = group_idx
                for i in range(len(blks)):
                    if adj[blk2int[curr_blk], i] == 1:
                        queue.append(int2blk[i])
            group_idx += 1

        self._name2alignment_group = name2alignment_group
        return name2alignment_group
    
    @property
    def name2alignment_group_color(self) -> dict[str, int]:
        if hasattr(self, '_name2alignment_group_color'):
            return self._name2alignment_group_color

        colors = ["red", "green", "yellow"]
        name2alignment_group_color = {}
        for blk in self.block_info:
            if blk.partner_indices.__len__() > 0:
                name2alignment_group_color[blk.name] = colors[self.name2alignment_group[blk.name] % len(colors)]
            else:
                name2alignment_group_color[blk.name] = "blue"

        self._name2alignment_group_color = name2alignment_group_color
        return name2alignment_group_color