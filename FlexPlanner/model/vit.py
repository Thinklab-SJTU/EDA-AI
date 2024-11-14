import torch
from torch import nn

from einops import rearrange, repeat
from einops.layers.torch import Rearrange

from typing import Tuple, Union
import numpy as np
from .transformer import Transformer as TransformerGNN
from tianshou.data import Batch



# helpers

def pair(t):
    return t if isinstance(t, tuple) else (t, t)

# classes

class FeedForward(nn.Module):
    def __init__(self, dim, hidden_dim, dropout = 0.):
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(dim),
            nn.Linear(dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, dim),
            nn.Dropout(dropout)
        )

    def forward(self, x):
        return self.net(x)

class Attention(nn.Module):
    def __init__(self, dim, heads = 8, dim_head = 64, dropout = 0.):
        super().__init__()
        inner_dim = dim_head *  heads
        project_out = not (heads == 1 and dim_head == dim)

        self.heads = heads
        self.scale = dim_head ** -0.5

        self.norm = nn.LayerNorm(dim)

        self.attend = nn.Softmax(dim = -1)
        self.dropout = nn.Dropout(dropout)

        self.to_qkv = nn.Linear(dim, inner_dim * 3, bias = False)

        self.to_out = nn.Sequential(
            nn.Linear(inner_dim, dim),
            nn.Dropout(dropout)
        ) if project_out else nn.Identity()

    def forward(self, x):
        x = self.norm(x)

        qkv = self.to_qkv(x).chunk(3, dim = -1)
        q, k, v = map(lambda t: rearrange(t, 'b n (h d) -> b h n d', h = self.heads), qkv)

        dots = torch.matmul(q, k.transpose(-1, -2)) * self.scale

        attn = self.attend(dots)
        attn = self.dropout(attn)

        out = torch.matmul(attn, v)
        out = rearrange(out, 'b h n d -> b n (h d)')
        return self.to_out(out)

class Transformer(nn.Module):
    def __init__(self, dim, depth, heads, dim_head, mlp_dim, dropout = 0.):
        super().__init__()
        self.norm = nn.LayerNorm(dim)
        self.layers = nn.ModuleList([])
        for _ in range(depth):
            self.layers.append(nn.ModuleList([
                Attention(dim, heads = heads, dim_head = dim_head, dropout = dropout),
                FeedForward(dim, mlp_dim, dropout = dropout)
            ]))

    def forward(self, x):
        for attn, ff in self.layers:
            x = attn(x) + x
            x = ff(x) + x

        return self.norm(x)

class ViT(nn.Module):
    """
    This is the shared encoder, implemented by Vision Transformer (ViT).
    """
    def __init__(self, input_channels:int, graph:bool, final_shape:tuple[int, int, int], 
                 image_size:Union[int, Tuple[int,int]], patch_size:Union[int, Tuple[int,int]], 
                 dim, depth, heads, mlp_dim, dim_head, dropout, emb_dropout):
        super().__init__()
        assert final_shape[1] == final_shape[2], f"final_shape[1] != final_shape[2], got {final_shape[1]} and {final_shape[2]}"

        image_height, image_width = pair(image_size)
        patch_height, patch_width = pair(patch_size)

        assert image_height % patch_height == 0 and image_width % patch_width == 0, 'Image dimensions must be divisible by the patch size.'

        self.input_channels = input_channels
        self.final_shape = final_shape
        self.last_channel = np.prod(final_shape)

        num_patches = (image_height // patch_height) * (image_width // patch_width)
        patch_dim = self.input_channels * patch_height * patch_width

        pool = 'mean'
        assert pool in {'cls', 'mean'}, 'pool type must be either cls (cls token) or mean (mean pooling)'

        self.to_patch_embedding = nn.Sequential(
            Rearrange('b c (h p1) (w p2) -> b (h w) (p1 p2 c)', p1 = patch_height, p2 = patch_width),
            nn.LayerNorm(patch_dim),
            nn.Linear(patch_dim, dim),
            nn.LayerNorm(dim),
        )

        self.pos_embedding = nn.Parameter(torch.randn(1, num_patches + 1, dim))
        self.dropout = nn.Dropout(emb_dropout)

        self.transformer = Transformer(dim, depth, heads, dim_head, mlp_dim, dropout)

        self.pool = pool
        self.to_latent = nn.Identity()

        self.mlp_head = nn.Linear(dim, self.last_channel)


        # graph
        self.graph = graph
        if graph:
            in_dim = 7
            hid_dim = 32
            out_dim = self.last_channel
            self.transformer_gnn = TransformerGNN(in_dim, hid_dim, out_dim, n_layers=2)

    def forward(self, stacked_mask:torch.Tensor, graph_data_batch:Batch) -> torch.Tensor:
        """
        @params:
            stacked_mask: shape = [B, C, H, W]
        @returns:
            shape = [B, C]
        """
        x = self.to_patch_embedding(stacked_mask)
        b, n, _ = x.shape
        

        x += self.pos_embedding[:, :n]

        x = self.dropout(x)

        x = self.transformer(x)

        x = x.mean(dim = 1) if self.pool == 'mean' else x[:, 0]

        x = self.to_latent(x)
        output = self.mlp_head(x)

        if self.graph:
            graph_data = []
            for k in graph_data_batch.keys():
                if k in ['x', 'y', 'z', 'w', 'h', 'area', 'placed']:
                    graph_data.append(graph_data_batch[k])
            graph_data = torch.stack(graph_data, dim=-1).to(output.device)

            adj_mat = graph_data_batch['adj_mat_mov'].to(output.device)
            curr_node = graph_data_batch['idx'].long()
            order = graph_data_batch['order'].long() if 'order' in graph_data_batch.keys() else None

            emb_global, emb_local = self.transformer_gnn(graph_data, adj_mat, curr_node, order)
            output = output + emb_global + emb_local

        return output

if __name__ == '__main__':
    input_channels = 6
    final_shape = (16, 8, 8)
    image_size = 64
    patch_size = 8
    dim = 32
    depth = 2
    heads = 4
    mlp_dim = 64
    dim_head = 32
    dropout = 0.0
    emb_dropout = 0.0

    batch_size = 4

    vit_cls = ViT
    if vit_cls is ViT:
        print('Using ViT')

    vit = ViT(input_channels, final_shape, image_size, patch_size, dim, depth, heads, mlp_dim, dim_head, dropout, emb_dropout)
    x = torch.rand((batch_size, input_channels, image_size, image_size))

    o = vit(x)
    print(o.shape)
    # torch.save(vit.state_dict(), 'vit.pt')
