import torch
from torch import nn
from einops import rearrange, repeat
from tianshou.data import Batch
import config

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
