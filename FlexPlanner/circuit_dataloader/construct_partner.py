from typing import List, Dict
from fp_env import Block, Net, FPInfo
from collections import OrderedDict
import numpy as np
import pandas as pd

def construct_partner_blk(fp_info:FPInfo, num_alignment:int, alignment_sort:str) -> pd.DataFrame:
    """
    Create alignment partners, return a dataframe with columns: blk0, blk1.
    """
    if num_alignment is None:
        num_alignment = 0
    
    block_info = fp_info.block_info
    df = pd.DataFrame(columns=["blk0", "blk1"])

    blk_idx_with_partner = set()

    if alignment_sort == "area":
        for _ in range(num_alignment):
            # get the block with the largest area and no partner
            tmp0 = [blk for blk in block_info if blk.idx not in blk_idx_with_partner]
            if len(tmp0) == 0:
                break
            tmp0.sort(key=lambda x: x.area, reverse=True)
            blk0 = tmp0[0]

            # get the block with the largest area and no partner and different layer
            tmp1 = [blk for blk in block_info if blk.idx not in blk_idx_with_partner and blk.grid_z != blk0.grid_z]
            if len(tmp1) == 0:
                break
            tmp1.sort(key=lambda x: x.area, reverse=True)
            blk1 = tmp1[0]

            if blk0.grid_z > blk1.grid_z:
                blk0, blk1 = blk1, blk0

            blk_idx_with_partner.add(blk0.idx)
            blk_idx_with_partner.add(blk1.idx)

            line = {
                "blk0": blk0.name,
                "blk1": blk1.name,
            }
            df = pd.concat([df, pd.DataFrame([line])])

    else:
        raise NotImplementedError(f"alignment_sort={alignment_sort} is not implemented.")
    
    return df