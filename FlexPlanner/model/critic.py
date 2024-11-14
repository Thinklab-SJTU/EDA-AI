import torch
import torchvision
from torch import nn
from tianshou.data import Batch, to_torch_as
from .shared_encoder import SharedEncoder
import numpy as np
from einops.layers.torch import Rearrange
import config
from .sequence_encoder import SequenceEncoderCNN as SequenceEncoderCNN2
from .layer_decider import DieEmbedding
from .transformer2 import Transformer2

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

class Critic(nn.Module):
    def __init__(self, 
            len_episode:int, shared_encoder:SharedEncoder, norm_wiremask:bool, 
            input_partner_die:bool, input_alignment_mask:bool, input_next_block:int, 
            num_die:int, input_sequence_critic:str, input_die_critic:bool, reduced_dim_critic:int, set_vision_to_zero:int, set_canvas_to_zero:int,
        ):
        super().__init__()
        self.register_buffer("redundancy", torch.tensor(0.0))
        self.norm_wiremask = norm_wiremask
        self.input_partner_die = input_partner_die
        self.input_alignment_mask = input_alignment_mask
        self.input_next_block = input_next_block
        self.input_sequence_critic = input_sequence_critic
        self.input_die_critic = input_die_critic
        self.reduced_dim_critic = reduced_dim_critic
        self.set_vision_to_zero = set_vision_to_zero
        self.set_canvas_to_zero = set_canvas_to_zero

        final_input_dim = 0

        # visual information backbone
        self.shared_encoder = shared_encoder
        self.visual_info_projector = nn.Linear(shared_encoder.last_channel, self.reduced_dim_critic) if self.reduced_dim_critic > 0 else nn.Identity()
        visual_embedding_dim = self.reduced_dim_critic if self.reduced_dim_critic > 0 else shared_encoder.last_channel
        final_input_dim += visual_embedding_dim


        # time embedding for current step
        time_embedding_dim = visual_embedding_dim
        self.time_embedding = nn.Embedding(len_episode+1, time_embedding_dim)
        

        # input sequence into critic
        if self.input_sequence_critic is not None:
            seq_in_dim = 8
            seq_hid_dim = 16
            seq_out_dim = 32
            final_input_dim += seq_out_dim*num_die
            if self.input_sequence_critic == "CNN":
                self.seq_enc = SequenceEncoderCNN(seq_in_dim, seq_hid_dim, seq_out_dim, num_die)
            elif self.input_sequence_critic == "Transformer2":
                self.seq_enc = Transformer2(seq_in_dim, seq_hid_dim, seq_out_dim, nhead=4, num_layers=2)
            elif self.input_sequence_critic == "CNN2":
                self.seq_enc = SequenceEncoderCNN2(seq_in_dim, seq_hid_dim, seq_out_dim, num_die)
            else:
                raise NotImplementedError(f"input_sequence_critic={self.input_sequence_critic} is not supported")


        # input die into critic
        if self.input_die_critic:
            final_input_dim += 2*num_die # two parts: 1) die embedding, 2) num_blk_without_placing_order
            self.die_embedding = nn.Embedding(num_die, num_die)
            # self.die_embedding = DieEmbedding(num_die, num_die, trainable=False)


        # final fusion module, output state value
        fusion_hid_dim = final_input_dim // 2 if self.reduced_dim_critic > 0 else self.shared_encoder.last_channel // 2 # be compatible to previous version
        self.fusion_module = nn.Sequential(
            nn.Linear(final_input_dim, fusion_hid_dim),
            nn.ReLU(),
            nn.Linear(fusion_hid_dim, 1)
        )

        
    @property
    def device(self) -> torch.device:
        return self.redundancy.device
    

    def forward(self, obs:Batch) -> torch.Tensor:
        """
        return the state value of the current state, shape is [B]
        """
        device = self.device
        step = torch.from_numpy(obs["step"]).to(device) # [B]
        canvas = obs["canvas"].to(device) # [B, 2, H, W]
        wiremask = obs["wiremask"].to(device) # [B, H, W]
        position_mask = obs["position_mask"].to(device) # [B, H, W]
        wiremask_next = obs["wiremask_next"].to(device) # [B, H, W]
        position_mask_next = obs["position_mask_next"].to(device) # [B, H, W]

        if self.set_canvas_to_zero:
            canvas = torch.zeros_like(canvas)

        if self.norm_wiremask:
            num_grid = wiremask.shape[-1]
            num_net = torch.from_numpy(obs["num_net"]).to(device=device, dtype=wiremask.dtype) # [B]
            wiremask = wiremask / num_net
            wiremask_next = wiremask_next / num_net
            wiremask_next_another = wiremask_next_another / num_net if wiremask_next_another is not None else None

        if self.input_partner_die:
            num_die = canvas.shape[1]
            stacked_mask = [canvas[:,d] for d in range(num_die)] + [wiremask, position_mask] # [B, D+2, H, W]
            stacked_mask = torch.stack(stacked_mask, dim=1)
        else:
            layer_idx = obs["layer_idx"] # [B]
            b = np.arange(canvas.shape[0]) # [B]
            stacked_mask = torch.stack([canvas[b,layer_idx], wiremask, position_mask], dim=1) # [B, 5, H, W]
        
        if self.input_next_block == 1:
            stacked_mask = torch.cat([stacked_mask, wiremask_next.unsqueeze(1), position_mask_next.unsqueeze(1)], dim=1)

        if self.input_alignment_mask:
            alignment_mask = obs["alignment_mask"].to(device) # [B, H, W]
            stacked_mask = torch.cat([stacked_mask, alignment_mask.unsqueeze(1)], dim=1)
            
            
        # set vision to zero, for ablation study
        if self.set_vision_to_zero:
            stacked_mask = torch.zeros_like(stacked_mask)

        # visual information backbone
        global_feat = self.shared_encoder(stacked_mask, to_torch_as(obs['graph_data'], stacked_mask) if self.shared_encoder.graph else None) # [B, C]
        global_feat = self.visual_info_projector(global_feat) # [B, C]

        # time embedding for current step
        time_feat = self.time_embedding(step) # [B, C]
        fusion_feat = global_feat + time_feat # add time embedding and global feature

        # input sequence into critic
        if self.input_sequence_critic is not None:
            seq_feat:Batch = obs["sequence_feature"]
            seq_feat = to_torch_as(seq_feat, fusion_feat)
            seq_feat = self.seq_enc(seq_feat)
            fusion_feat = torch.cat([fusion_feat, seq_feat], dim=-1)
        
        # input die into critic
        if self.input_die_critic:
            layer_idx = torch.from_numpy(obs["layer_idx"]).to(device=device, dtype=torch.long) # [B]
            num_blk_without_placing_order = torch.from_numpy(obs["num_blk_without_placing_order"]).to(device=device, dtype=stacked_mask.dtype) # [B, num_die]
            die_embedding = self.die_embedding(layer_idx)
            die_feat = torch.cat([die_embedding, num_blk_without_placing_order], dim=-1)
            fusion_feat = torch.cat([fusion_feat, die_feat], dim=-1)

        state_value = self.fusion_module(fusion_feat).squeeze(-1) # [B]
        return state_value
