import math
import numpy as np
from .terminal import Terminal

class Block:
    def __init__(self, x:float, y:float, z:int, w:float, h:float, 
                 name:str, preplaced:bool, virtual:bool, power:float,
                 ):
        """
        x: bottom-left x
        y: bottom-left y
        z: die/layer index
        """
        # basic properties
        self.x, self.y, self.w, self.h = x, y, w, h
        self.z = self.grid_z = z
        self.name = name
        self.preplaced = bool(preplaced)
        self.virtual = bool(virtual)
        self.power = power

        # alignment, use full idx
        self.partner_indices:list[int] = []
        self.alignment_areas:dict[int,int] = dict()
        self.original_alignment_areas:dict[int,float] = dict()

        # adjacent terminals
        self.adjacent_terminals:list[Terminal] = []

        # adjacent blocks, use full idx
        self.adjacent_blocks:list[int] = []

        # init properties
        self.init_w, self.init_h = w, h


        # connected_nets
        self.connected_nets = []

        # reset
        self.reset()
    

    def set_idx(self, full_idx:int, self_part_idx:int):
        """self_part_idx is either preplaced_idx or movable_idx."""
        self.idx = full_idx
        if self.preplaced:
            self.preplaced_idx = self_part_idx
        else:
            self.movable_idx = self_part_idx


    def set_grid_wh(self, grid_width:float, grid_height:float):
        """set w and h based on grid_width and grid_height."""
        self.grid_w = round(max(1, self.w / grid_width)) # at least 1 grid
        self.grid_h = round(max(1, self.h / grid_height)) # at least 1 grid
        self.init_grid_w, self.init_grid_h = self.grid_w, self.grid_h
    

    def set_grid_xy(self, grid_width:float, grid_height:float, x_grid_num:int, y_grid_num:int):
        """set grid_x and grid_y based on x, y, grid_width and grid_height."""
        self.grid_x = round(min(x_grid_num - self.grid_w, self.x / grid_width)) # at most x_grid_num - block.grid_w
        self.grid_y = round(min(y_grid_num - self.grid_h, self.y / grid_height)) # at most y_grid_num - block.grid_h
    

    def reset(self):
        """reset placed to preplaced."""
        self.placed = self.preplaced
        self.w, self.h = self.init_w, self.init_h
        if hasattr(self, "init_grid_w") and hasattr(self, "init_grid_h"):
            self.grid_w, self.grid_h = self.init_grid_w, self.init_grid_h

    
    def set_ratio(self, ratio:float, x_grid_num:int, y_grid_num:int):
        """ratio is w/h."""
        # ratio = grid_w/grid_h, reset grid_w and grid_h based on ratio, keep grid_area unchanged
        original_grid_area = self.grid_area
        self.grid_w = max(1, round(math.sqrt(original_grid_area * ratio)))
        self.grid_h = max(1, round(self.grid_w / ratio))

        # avoid out of boundary
        if self.grid_w > x_grid_num:
            self.grid_w = x_grid_num
            self.grid_h = max(1, round(original_grid_area / self.grid_w))
            ratio = self.grid_w / self.grid_h

        if self.grid_h > y_grid_num:
            self.grid_h = y_grid_num
            self.grid_w = max(1, round(original_grid_area / self.grid_h))
            ratio = self.grid_w / self.grid_h

        # ratio = w/h, reset w and h based on ratio, keep area unchanged
        self.w = math.sqrt(self.area * ratio)
        self.h = self.w / ratio

        assert self.grid_w <= x_grid_num and self.grid_h <= y_grid_num, "[ERROR] block {} is too large: grid_w={}, grid_h={}, x_grid_num={}, y_grid_num={}".format(self.name, self.grid_w, self.grid_h, x_grid_num, y_grid_num)
        

    def rotate(self):
        """rotate the block."""
        self.w, self.h = self.h, self.w
        self.grid_w, self.grid_h = self.grid_h, self.grid_w
    

    def place(self, grid_x:int, grid_y:int):
        """place the block to grid_x, grid_y."""
        self.placed = True
        self.grid_x, self.grid_y = grid_x, grid_y


    @property
    def area(self) -> float:
        """return original area, not grid area"""
        return self.w * self.h
    

    @property
    def grid_area(self) -> int:
        """return grid area, not original area"""
        return self.grid_w * self.grid_h


    def __repr__(self) -> str:
        self_part_idx_str = "preplaced_idx={}".format(self.preplaced_idx) if self.preplaced else "movable_idx={}".format(self.movable_idx)
        return "Block(name={}, id={}, {}, placed={}, grid_x={}, grid_y={}, grid_z={}, grid_w={}, grid_h={}, grid_area={}, alignment_partners={})".format(
            self.name,
            self.idx, self_part_idx_str,
            self.placed,
            self.grid_x, self.grid_y, self.grid_z, self.grid_w, self.grid_h, self.grid_area, self.partner_indices,
        )


    def set_partner(self, partner_idx:int, alignment_area:int, original_alignment_area:float):
        """
        set partner block idx, full idx.
        alignment_area: alignment area between two blocks in discrete grid.
        """
        if partner_idx not in self.partner_indices:
            self.partner_indices.append(partner_idx)
            self.alignment_areas[partner_idx] = alignment_area
            self.original_alignment_areas[partner_idx] = original_alignment_area
        return
    

    def set_adjacent_terminal(self, terminal:Terminal):
        """set adjacent terminal."""
        if terminal not in self.adjacent_terminals:
            self.adjacent_terminals.append(terminal)
        return
    

    def calc_distance_adjacent_terminal(self, tml:Terminal) -> float:
        """
        Return the distance to the terminal.
        This distance is not normalized.
        The lower, the better.
        """
        assert self.placed, "[ERROR] block {} is not placed.".format(self.name)
        x,y = tml.grid_x, tml.grid_y
        # edge left
        x1 = self.grid_x
        y1,y2 = self.grid_y, self.grid_y + self.grid_h
        x_dis = abs(x1 - x)
        if y < y1:
            y_dis = y1 - y
        elif y > y2:
            y_dis = y - y2
        else:
            y_dis = 0
        dis_left = min(x_dis, y_dis)

        # edge right
        x1 = self.grid_x + self.grid_w
        y1,y2 = self.grid_y, self.grid_y + self.grid_h
        x_dis = abs(x1 - x)
        if y < y1:
            y_dis = y1 - y
        elif y > y2:
            y_dis = y - y2
        else:
            y_dis = 0
        dis_right = min(x_dis, y_dis)

        # edge top
        x1,x2 = self.grid_x, self.grid_x + self.grid_w
        y1 = self.grid_y
        y_dis = abs(y1 - y)
        if x < x1:
            x_dis = x1 - x
        elif x > x2:
            x_dis = x - x2
        else:
            x_dis = 0
        dis_top = min(x_dis, y_dis)

        # edge bottom
        x1,x2 = self.grid_x, self.grid_x + self.grid_w
        y1 = self.grid_y + self.grid_h
        y_dis = abs(y1 - y)
        if x < x1:
            x_dis = x1 - x
        elif x > x2:
            x_dis = x - x2
        else:
            x_dis = 0
        dis_bottom = min(x_dis, y_dis)

        dis_min = min(dis_left, dis_right, dis_top, dis_bottom)
        return dis_min


    def set_adjacent_block(self, block_idx:int):
        """set adjacent block."""
        if block_idx not in self.adjacent_blocks:
            self.adjacent_blocks.append(block_idx)
        return
    

    def calc_reward_adjacent_block(self, block, reward_class_adjacent_block:int, num_x_grid:int, num_y_grid:int) -> tuple[float,float]:
        """
        Return the adjacent length and reward.
        """
        assert self.placed, "[ERROR] block {} is not placed.".format(self.name)
        assert block.placed, "[ERROR] block {} is not placed.".format(block.name)

        def calc_overlap_1d(x1:int, x2:int, w1:int, w2:int) -> int:
            """calculate the overlap area for 1D."""
            left, right = max(x1, x2), min(x1 + w1, x2 + w2)
            return max(0, right - left)

        x1,y1,w1,h1 = self.grid_x, self.grid_y, self.grid_w, self.grid_h
        x2,y2,w2,h2 = block.grid_x, block.grid_y, block.grid_w, block.grid_h
        flag = True

        if x2 + w2 == x1:
            adj_len = calc_overlap_1d(y1, y2, h1, h2)
        elif x1 + w1 == x2:
            adj_len = calc_overlap_1d(y1, y2, h1, h2)
        elif y2 + h2 == y1:
            adj_len = calc_overlap_1d(x1, x2, w1, w2)
        elif y1 + h1 == y2:
            adj_len = calc_overlap_1d(x1, x2, w1, w2)
        else:
            adj_len = 0
            flag = False # adjacent length is 0


        # calculation for reward
        # blk is adjacent to blk
        if flag:
            if reward_class_adjacent_block in [0,1,2]:
                reward = adj_len
            else:
                raise NotImplementedError("[ERROR] reward_class_adjacent_block={} is not implemented.".format(reward_class_adjacent_block))
        
        # blk is not adjacent to blk
        else:
            x_overlap, y_overlap = calc_overlap_1d(x1, x2, w1, w2), calc_overlap_1d(y1, y2, h1, h2)

            # block overlap
            if x_overlap > 0 and y_overlap > 0:
                if reward_class_adjacent_block == 0:
                    reward = 0
                elif reward_class_adjacent_block in [1,2]:
                    overlap = x_overlap * y_overlap
                    reward = -(overlap**0.5) - (num_x_grid + num_y_grid)/2 # minus a const
                else:
                    raise NotImplementedError("[ERROR] reward_class_adjacent_block={} is not implemented.".format(reward_class_adjacent_block))
            
            # one side overlap, but not adjacent to blk
            elif x_overlap > 0:
                y_dis = min( abs(y1-y2), abs(y1+h1-y2-h2), abs(y1-y2-h2), abs(y1+h1-y2) )
                if reward_class_adjacent_block == 0:
                    reward = 0
                elif reward_class_adjacent_block == 1:
                    reward = -y_dis
                elif reward_class_adjacent_block == 2:
                    reward = (1 - y_dis/num_y_grid) * x_overlap
                else:
                    raise NotImplementedError("[ERROR] reward_class_adjacent_block={} is not implemented.".format(reward_class_adjacent_block))
                

            elif y_overlap > 0:
                x_dis = min( abs(x1-x2), abs(x1+w1-x2-w2), abs(x1-x2-w2), abs(x1+w1-x2) )
                if reward_class_adjacent_block == 0:
                    reward = 0
                elif reward_class_adjacent_block == 1:
                    reward = -x_dis
                elif reward_class_adjacent_block == 2:
                    reward = (1 - x_dis/num_x_grid) * y_overlap
                else:
                    raise NotImplementedError("[ERROR] reward_class_adjacent_block={} is not implemented.".format(reward_class_adjacent_block))

            # not close to blk, and neither 1d nor 2d overlap
            else:
                if reward_class_adjacent_block == 0:
                    reward = 0
                elif reward_class_adjacent_block in [1,2]:
                    x_dis = min( abs(x1-x2), abs(x1+w1-x2-w2), abs(x1-x2-w2), abs(x1+w1-x2) )
                    y_dis = min( abs(y1-y2), abs(y1+h1-y2-h2), abs(y1-y2-h2), abs(y1+h1-y2) )
                    distance = min(x_dis, y_dis)
                    reward = -distance
                else:
                    raise NotImplementedError("[ERROR] reward_class_adjacent_block={} is not implemented.".format(reward_class_adjacent_block))

        return adj_len, reward     