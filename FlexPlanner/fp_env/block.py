import math
import numpy as np
from .terminal import Terminal

class Block:
    def __init__(self, x:float, y:float, z:int, w:float, h:float, 
                 name:str, preplaced:bool, virtual:bool,
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

        # alignment, use full idx
        self.partner_indices:list[int] = []
        self.alignment_areas:dict[int,int] = dict()
        self.original_alignment_areas:dict[int,float] = dict()


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
