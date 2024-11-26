import torch
from torch import nn
from tianshou.data import Batch
from einops.layers.torch import Rearrange, Reduce
import config
from einops import rearrange, repeat


class DieEmbedding(nn.Module):
    def __init__(self, num_die:int, dim:int, trainable:bool=True):
        super().__init__()
        self.trainable = trainable
        if self.trainable:
            self.die_embedding = nn.Embedding(num_die+1, dim)
        else:
            t = torch.arange(num_die+1).float()
            t[0] = -1.0
            t = repeat(t, 'n -> n c', c=dim)
            self.register_buffer('die_embedding', t)
    
    def forward(self, die_idx:torch.LongTensor) -> torch.Tensor:
        """
        @param die_idx: [B], B is batch_size
        @return: [B, C], C is dim
        """
        if self.trainable:
            return self.die_embedding(die_idx)
        else:
            return self.die_embedding[die_idx]
        

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
        """in_channels and out_channels are dim for each die, not the total."""
        super().__init__(in_channels, num_die)
        self.rearrange = Rearrange('B L N C -> B (L C) N')
        self.seq_encoder = nn.Sequential(
            nn.Conv1d(num_die*in_channels, num_die*hidden_channels, 7, 1, 3, groups=num_die),
            nn.ReLU(),
            nn.Conv1d(num_die*hidden_channels, num_die*hidden_channels, 5, 1, 2, groups=num_die),
            nn.ReLU(),
            nn.Conv1d(num_die*hidden_channels, num_die*out_channels, 3, 1, 1, groups=num_die),
            nn.AdaptiveMaxPool1d(1), # [B, num_die*out_channels, 1]
            Rearrange('b c 1 -> b c'), # [B, num_die*out_channels]
        )


class SequenceEncoderCNN2(SequenceEncoderBase):
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



class Transformer2(nn.Module):
    """
    This is the backbone to process the sequence feature, i.e., the blocks placing orders in multiple dies.
    """
    def __init__(self, in_dim:int, hid_dim:int, out_dim:int, nhead:int=4, num_layers:int=2):
        super().__init__()
        self.nhead = nhead

        self.fc1 = nn.Linear(in_dim, hid_dim)
        self.fc2 = nn.Linear(hid_dim, out_dim)

        enc_layer = nn.TransformerEncoderLayer(hid_dim, nhead, hid_dim*2, dropout=0.0, batch_first=True)
        dec_layer = nn.TransformerDecoderLayer(hid_dim, nhead, hid_dim*2, dropout=0.0, batch_first=True)

        self.encoder = nn.TransformerEncoder(enc_layer, num_layers)
        self.decoder = nn.TransformerDecoder(dec_layer, num_layers)



    # def forward(self, seq:torch.Tensor, curr_blk:torch.Tensor, place_order:torch.Tensor) -> torch.Tensor:
    def forward(self, x:Batch) -> torch.Tensor:
        """
        @params:
            seq: [B, num_die, num_blk, in_channels]
            curr_blk: [B, in_channels]
            place_order: [B, num_die, num_blk]
        @return:
            output: [B, num_die*out_dim]
        """
        # get data
        seq = torch.stack([x[key] for key in config.sequence_feature_keys], dim=-1) # [B, num_die, num_blk, in_channels]
        curr_blk = torch.stack([x["last_placed_block"][key] for key in config.sequence_feature_keys], dim=-1) # [B, in_channels]
        place_order = x.sequence.long() # [B, num_die, num_blk]

        # preprocess
        seq = self.fc1(seq)
        curr_blk = self.fc1(curr_blk)

        # mask
        mask = torch.where(place_order == -1, 1, 0).to(seq.device) # [B N L]
        m1 = rearrange(mask, 'B N L -> B N L 1')
        m2 = rearrange(mask, 'B N L -> B N 1 L')
        mask = m1 + m2
        mask[:,:,torch.arange(mask.size(2)),torch.arange(mask.size(3))] = 0
        mask = mask.bool() # [B N L L]
        mask = rearrange(mask, 'B N L1 L2 -> (B N) L1 L2')
        mask = repeat(mask, 'B L1 L2 -> B H L1 L2', H=self.nhead)
        mask = rearrange(mask, 'B H L1 L2 -> (B H) L1 L2')

        # encoder
        num_die = seq.size(1)
        src = rearrange(seq, 'B N L C -> (B N) L C')
        memory = self.encoder(src, mask)

        # decoder
        tgt = repeat(curr_blk, 'B C -> B N 1 C', N=num_die)
        tgt = rearrange(tgt, 'B N 1 C -> (B N) 1 C')
        mask = torch.where(place_order == -1, 1, 0).to(tgt.device) # [B N L]
        mask = rearrange(mask, 'B N L -> (B N) 1 L')
        mask = repeat(mask, 'B 1 L -> B H 1 L', H=self.nhead)
        mask = rearrange(mask, 'B H 1 L -> (B H) 1 L')
        mask = mask.bool()
        output = self.decoder(tgt, memory, memory_mask=mask)
        output = rearrange(output, '(B N) 1 C -> B N C', N=num_die)

        # output
        output = self.fc2(output)
        output = rearrange(output, 'B N C -> B (N C)')

        return output