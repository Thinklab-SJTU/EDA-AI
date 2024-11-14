import torch
from torch import nn
from einops import rearrange, repeat


class PositionalEncoding(nn.Module):
    def __init__(self, max_len:int, dim:int, trainable:bool=True):
        super().__init__()
        if trainable:
            self.pos_enc = nn.Parameter(torch.randn(max_len, dim))
        else:
            pos_enc = torch.zeros(max_len, dim)
            pos = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
            div_term = torch.exp(torch.arange(0, dim, 2).float() * (-torch.log(torch.tensor(10000.0)) / dim))
            pos_enc[:, 0::2] = torch.sin(pos * div_term)
            pos_enc[:, 1::2] = torch.cos(pos * div_term)
            self.register_buffer('pos_enc', pos_enc)

    def forward(self, x:torch.Tensor, order:torch.LongTensor) -> torch.Tensor:
        """order: [B,L]"""
        if order is None:
            return x
        
        batch_size, seq_len = order.shape
        pe = self.pos_enc[order.flatten()]
        pe = pe.view(batch_size, seq_len, -1)
        return pe + x

class Transformer(nn.Module):
    """
    This is the backbone to process graph feature in the shared encoder.
    """
    def __init__(self, in_dim:int, hid_dim:int, out_dim:int, n_layers:int):
        super().__init__()
        self.fc0 = nn.Linear(in_dim, hid_dim)
        self.fc1 = nn.Linear(hid_dim, out_dim)
        self.fc2 = nn.Linear(hid_dim, out_dim)
        self.nhead = 4

        enc_layer = nn.TransformerEncoderLayer(d_model=hid_dim, nhead=self.nhead, dim_feedforward=hid_dim*2, dropout=0.0, batch_first=True)
        self.transformer = nn.TransformerEncoder(enc_layer, num_layers=n_layers)
        self.invalid_node = nn.Parameter(torch.randn(1, hid_dim))

        self.pe = PositionalEncoding(300, hid_dim, trainable=False)

    def forward(self, x:torch.Tensor, adj_mat:torch.IntTensor, idx:torch.LongTensor, order:torch.LongTensor) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Input:
            x: the input tensor, is a torch.Tensor, shape = [B, L, D].
        Output:
            output: global and local embeddings, both are torch.Tensor, shape = [B, D].
        """
        x = self.fc0(x)
        x = self.pe(x, order)
        mask = (1-adj_mat).bool()
        mask = repeat(mask, 'b i j -> b n i j', n=self.nhead)
        mask = rearrange(mask, 'b n i j -> (b n) i j')
        x = self.transformer(x, mask)

        if (idx == -1).sum() > 0:
            pass
        
        emb_local = x[torch.arange(x.size(0)), idx]
        emb_local[idx == -1] = self.invalid_node
        emb_global = x.mean(dim=1)

        emb_local = self.fc1(emb_local)
        emb_global = self.fc2(emb_global)
        return emb_global, emb_local
