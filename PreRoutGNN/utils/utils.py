import random 
import json
import numpy as np
import random
import os
import torch
import pickle
import time
import argparse
import shutil
import dgl
import copy
import functools

def record_time(func):
    @functools.wraps(func)
    def run(*args, **kwds):
        torch.cuda.synchronize()
        s = time.time()
        ans = func(*args, **kwds)
        torch.cuda.synchronize()
        e = time.time()
        print(f"Running time for {func.__name__} = {e-s} (s)")
        return ans
    return run


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
    os.environ['PYTHONHASHSEED'] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    dgl.seed(seed)


def get_datetime():
    t = time.strftime('%Y-%m-%d-%H-%M-%S', time.localtime())
    return t


class Logger():
    def __init__(self,log_file_path) -> None:
        self.path = log_file_path
        with open(self.path,'w') as f:
            f.write(get_datetime() + "\n")
            print(get_datetime())
        return
    
    def log(self,content):
        content = str(content)
        with open(self.path,'a') as f:
            f.write(content + "\n")
            print(content)
        return


def mkdir(dir, rm=False):
    if os.path.isdir(dir) and rm:
        shutil.rmtree(dir)
    os.makedirs(dir)


def convert_dict_to_args(d):
    parser = argparse.ArgumentParser()
    for k,v in d.items():
        parser.add_argument(f'--{k}', default=v)
    return parser.parse_args()


def check_topo(g:dgl.heterograph):
    n_src, n_dst = g.edges(etype='net_out')
    c_src, c_dst = g.edges(etype='cell_out')
    g_homo = dgl.graph((
        torch.cat([n_src, c_src]),
        torch.cat([n_dst, c_dst]),
    ))
    try:
        topo_levels = dgl.topological_nodes_generator(g_homo)
        num_topo_levels = len(topo_levels)
        if num_topo_levels % 2 == 0:
            return True
        else:
            print(f"Number of topo levels is {num_topo_levels}, which must be even")
            return False
    except dgl.DGLError as e:
        # print(repr(e))
        print("Loop detected")
        return False


def get_from_dicts(key, *dicts):
    for d in dicts:
        if key in d.keys():
            return d[key]
    raise KeyError(key)


def merge_dicts(*dicts):
    res = {}
    for d in dicts:
        for k in d.keys():
            if k in res.keys():
                continue
            else:
                res[k] = d[k]
    return copy.deepcopy(res)


def write_model_structure(model, path:str):
    with open(path, 'w') as f:
        f.write(str(model))


def mkdir(dir, rm=False):
    if os.path.isdir(dir):
        if rm:
            shutil.rmtree(dir)
            os.makedirs(dir)
        else:
            pass
    else:
        os.makedirs(dir)

# @record_time
def to_device(g, ts, device):
    g = g.to(device)
    for k in ts.keys():
        if torch.is_tensor(ts[k]):
            ts[k] = ts[k].to(device)
        elif isinstance(ts[k], dgl.DGLGraph):
            ts[k] = ts[k].to(device)
        elif isinstance(ts[k], list):
            ts[k] = [i.to(device) for i in ts[k]]
    return g, ts


if __name__ == "__main__":
    print(get_datetime())