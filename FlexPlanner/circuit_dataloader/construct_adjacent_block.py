import fp_env
import pandas as pd

def construct_block_adjacent_block(fp_info:fp_env.FPInfo) -> pd.DataFrame:
    """
    return a dataframe with two columns: blk0, blk1
    """
    df = pd.DataFrame(columns=['blk0', 'blk1'])
    return df