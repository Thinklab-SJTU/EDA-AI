import torch
import pandas as pd

def write_floorplan(blk_xy:torch.Tensor, blk_wh:torch.Tensor, save_path:str):
    # blk_xy: (N, 2)
    # blk_wh: (N, 2)
    # save_path: str
    # return: None
    # write to save_path, with following format:
    # num_blk
    # x1 y1 w1 h1
    # x2 y2 w2 h2
    # ...
    # xn yn wn hn
    if not save_path.endswith('.csv'):
        num_blk = blk_xy.shape[0]
        with open(save_path, 'w') as f:
            f.write(str(num_blk) + '\n')
            for i in range(num_blk):
                x, y = blk_xy[i]
                w, h = blk_wh[i]
                f.write(f"{x} {y} {w} {h}\n")
    else:
        # save blk_xy and blk_wh to csv
        data = torch.cat([blk_xy, blk_wh], dim=1).detach().cpu().numpy()
        df = pd.DataFrame(data, columns=['x', 'y', 'w', 'h'])
        df.to_csv(save_path, index=False)
        
