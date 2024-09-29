import torch
from einops import rearrange, repeat
import math

@torch.no_grad()
def get_positional_encoding(x, num_freq_compents, max_L:float=1, frequency:torch.FloatTensor=None):
    """x.shape = [batch_size,]"""
    assert len(x.shape) == 1, f'len(x.shape) should be 1, but got x\'s shape is {x.shape}'
    assert frequency is None, f"frequency in get_positional_encoding should be None, but got {frequency}"
    
    L = num_freq_compents
    x_sin = repeat(x, "B -> B L", L=L).clone()
    x_cos = repeat(x, "B -> B L", L=L).clone()
    L = torch.full(size=(1,L), fill_value=2).pow(torch.arange(L).view(1,L)).to(dtype=x.dtype, device=x.device)
    x_sin = torch.sin(x_sin * math.pi * L / max_L)
    x_cos = torch.cos(x_cos * math.pi * L / max_L)
    x_embedding = torch.cat([x_sin, x_cos], dim=1)
    return x_embedding


@torch.no_grad()
def get_gaussian_encoding(x, num_freq_compents, max_L:float=1, frequency:torch.FloatTensor=None):
    """x.shape = [batch_size,]"""
    assert frequency is not None, f'frequency should not be None'
    assert len(frequency.shape) == 1, f"len(frequency.shape) should be 1, but got frequency's shape is {frequency.shape}"
    assert frequency.shape[0] == num_freq_compents, f"frequency's shape is {frequency.shape}, but num_freq_compents is {num_freq_compents}"
    assert len(x.shape) == 1, f"len(x.shape) should be 1, but got x's shape is {x.shape}"
    
    frequency = frequency.to(x.device).reshape(1, -1)
    x_sin = repeat(x, "B -> B L", L=num_freq_compents).clone()
    x_cos = repeat(x, "B -> B L", L=num_freq_compents).clone()
    x_sin = torch.sin(x_sin * 2 * math.pi * frequency / max_L)
    x_cos = torch.cos(x_cos * 2 * math.pi * frequency / max_L)
    x_embedding = torch.cat([x_sin, x_cos], dim=1)
    return x_embedding