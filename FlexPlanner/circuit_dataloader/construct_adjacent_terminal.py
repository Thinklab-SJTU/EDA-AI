import fp_env
import pandas as pd

def construct_block_adjacent_terminal(fp_info:fp_env.FPInfo) -> pd.DataFrame:
    """
    return a DataFrame with columns ['blk', 'tml']
    """
    df = pd.DataFrame(columns=['blk', 'tml'])
    return df
