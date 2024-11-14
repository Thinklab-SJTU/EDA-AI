import random 
import json
import numpy as np
import random
import os
import torch
import pickle
import time
import shutil
import functools
from typing import Any

tqdm_config = {
    "dynamic_ncols": True,
    "ascii": True,
}

def load_json(path):
    with open(path,'r') as f:
        res = json.load(f)
    return res


def save_json(obj, path:str):
    with open(path, 'w', encoding='utf8') as f:
        json.dump(obj, f, indent=4)


def load_pkl(path):
    with open(path, 'rb') as f:
        res = pickle.load(f)
        return res


def save_pkl(obj, path):
    with open(path, 'wb') as f:
        pickle.dump(obj, f)


def setup_seed(seed = 3407):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.cuda.manual_seed(seed)
    # https://zhuanlan.zhihu.com/p/73711222
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    os.environ["PYTHONHASHSEED"] = str(seed)


def get_datetime():
    t = time.strftime('%Y-%m-%d-%H-%M-%S', time.localtime())
    return t


def mkdir(dir:str, rm:bool=False):
    if os.path.isdir(dir):
        if rm:
            shutil.rmtree(dir)
            os.makedirs(dir)
        else:
            pass
    else:
        os.makedirs(dir)


def record_time(func):
    @functools.wraps(func)
    def run(*args, **kwds):
        torch.cuda.synchronize()
        s = time.time()
        ans = func(*args, **kwds)
        torch.cuda.synchronize()
        e = time.time()
        print("[INFO] Running time for [{:20}] = {:.6f} (s)".format(func.__name__, e-s))
        return ans
    return run



class MaskGenerator(torch.nn.Module):
    def __init__(self, mask_ratio:float, masked_padding_value:float) -> None:
        """
        @params:
            mask_ratio: ratio to masked
            masked_padding_value, replace original value with this
        """
        super().__init__()
        self.mask_ratio = mask_ratio
        self.masked_padding_value = masked_padding_value
    
    
    def forward(self, x:torch.Tensor):
        return self.generate_mask(x)

    
    def generate_mask(self, x:torch.Tensor):
        """
        @params:
            x: shape = [batch_size, seq_len, channel]
        @returns:
            masked_x
            mask: 1 means keep, 0 means masked
            mask_idx
            unmask_idx
        """
        mask_len = int(self.mask_ratio * x.shape[1])

        idx_shuffle = torch.rand_like(x).argsort(dim=1)

        # .sort() to keep their original relative order
        mask_idx = idx_shuffle[:, 0:mask_len, :].sort(dim=1)[0]
        unmask_idx = idx_shuffle[:, mask_len:, :].sort(dim=1)[0]
        restore_idx = torch.cat([mask_idx, unmask_idx], dim=1).argsort(dim=1)

        masked_part = torch.full_like(mask_idx, self.masked_padding_value, dtype=x.dtype)
        unmasked_part = x.gather(dim=1, index=unmask_idx)
        masked_x = torch.cat([masked_part, unmasked_part], dim=1).gather(dim=1, index=restore_idx)

        mask = torch.cat([torch.zeros_like(mask_idx), torch.ones_like(unmask_idx)], dim=1)
        mask = mask.gather(dim=1, index=restore_idx)

        return masked_x, mask, mask_idx, unmask_idx
    

class DummyTqdm:
    """A dummy tqdm class that keeps stats but without progress bar.

    It supports ``__enter__`` and ``__exit__``, update and a dummy
    ``set_postfix``, which is the interface that trainers use.

    .. note::

        Using ``disable=True`` in tqdm config results in infinite loop, thus
        this class is created. See the discussion at #641 for details.
    """

    def __init__(self, total: int, **kwargs: Any):
        self.total = total
        self.n = 0

    def set_postfix(self, **kwargs: Any) -> None:
        pass

    def update(self, n: int = 1) -> None:
        self.n += n

    def __enter__(self) -> "DummyTqdm":
        return self

    def __exit__(self, *args: Any, **kwargs: Any) -> None:
        pass


def is_power_of_2(n:int) -> bool:
    return (n & (n-1) == 0) and n != 0


def load_checkpoint_mismatch(model:torch.nn.Module, checkpoint:dict, allow_mismatch:bool=True) -> None:
    new_dict = model.state_dict()
    old_dict = checkpoint
    
    if allow_mismatch:
        for old_key in old_dict.keys():
            if old_key not in new_dict:
                print("[WARNING] old_key {} not in new_dict".format(old_key))
            elif old_dict[old_key].shape != new_dict[old_key].shape:
                print("[WARNING] key {} shape mismatch, old shape = {}, new shape = {}".format(old_key, old_dict[old_key].shape, new_dict[old_key].shape))  
            else:  
                new_dict[old_key] = old_dict[old_key]
        
        for new_key in set(new_dict.keys()) - set(old_dict.keys()):
            print("[WARNING] new_key {} not in old_dict".format(new_key))
                
    model.load_state_dict(new_dict)
    

def set_grad(model:torch.nn.Module, grad:bool) -> None:
    if model is None:
        return
    for param in model.parameters():
        param.requires_grad = grad

def set_grad_none(model:torch.nn.Module) -> None:
    for p in model.parameters():
        p.grad = None
