import torch
from torch import nn
from tianshou.data import Batch
from einops.layers.torch import Rearrange, Reduce
import config

class SequenceEncoderBase(nn.Module):
    def __init__(self, in_channels:int, num_die:int):
        super().__init__()
        self.rearrange = None
        self.seq_encoder = None
    
    def forward(self, seq_feat:Batch) -> torch.Tensor:
        """
        @param seq_feat: [B, L, N, C], L is num_die, N is num_blk
        @return: [B, L*out_channels], L is num_die
        """
        seq_feat = torch.stack([seq_feat[key] for key in config.sequence_feature_keys], dim=-1) # [B, L, N, C]
        x = self.rearrange(seq_feat)
        x = self.seq_encoder(x)
        return x


class SequenceEncoderCNN(SequenceEncoderBase):
    def __init__(self, in_channels:int, hidden_channels:int, out_channels:int, num_die:int):
        """
        This is the backbone to process the sequence feature, i.e., the blocks placing orders in multiple dies.
        """
        super().__init__(in_channels, num_die)
        self.rearrange = Rearrange('B L N C -> B (L C) N') # N is the number of blocks
        self.seq_encoder = nn.Sequential(
            nn.Conv1d(num_die*in_channels, num_die*hidden_channels, 3, 1, 1, groups=num_die),
            nn.ReLU(),
            nn.MaxPool1d(2),

            nn.Conv1d(num_die*hidden_channels, num_die*hidden_channels, 3, 1, 1, groups=num_die),
            nn.ReLU(),
            nn.MaxPool1d(2),

            nn.Conv1d(num_die*hidden_channels, num_die*out_channels, 3, 1, 1, groups=num_die),
            nn.AdaptiveAvgPool1d(1), # [B, num_die*out_channels, 1]
            Rearrange('b c 1 -> b c'), # [B, num_die*out_channels]
        )


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

    def forward(self, x:torch.Tensor) -> torch.Tensor:
        seq_len = x.size(1)
        pos_enc = self.pos_enc[:seq_len]
        pos_enc = pos_enc.unsqueeze(0).expand_as(x)
        return x + pos_enc
