import torch
import torchvision
from torch import nn
from einops import rearrange, repeat
from typing import Tuple
from tianshou.data import Batch, to_torch_as
from typing import Union
import numpy as np
from .shared_encoder import SharedEncoder
from .generator import InfoGANGenerator, TransposedConv
from .ratio_decider import RatioDecider
from .layer_decider import LayerDecider
from torch.distributions import constraints


def split(n:int, num:int) -> list[int]:
    """split n into num parts, each part should be bigger or equal to 1, and the sum of parts should be equal to n"""
    assert n >= num, f"n = {n}, num = {num}"
    parts = [n // num] * num
    for i in range(n % num):
        parts[i] += 1
    return parts


class LocalEncoder(nn.Module):
    def __init__(self, input_channel:int, output_channel:int):
        """keep image size unchanged (kernel_size=3, stride=1, padding=1)"""
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(input_channel, 4, 3, 1, 1),
            nn.ReLU(),
            nn.Conv2d(4, 8, 3, 1, 1),
            nn.ReLU(),
            nn.Conv2d(8, output_channel, 3, 1, 1),
        )

    
    def forward(self, stacked_mask:torch.Tensor) -> torch.Tensor:
        """stacked_mask.shape = [B, C, H, W]"""
        return self.conv(stacked_mask)
    

class Actor(nn.Module):
    def __init__(self, 
            num_act:int, shared_encoder:SharedEncoder, local_encoder:LocalEncoder, actor_update_shared_encoder:bool, ratio_decider:RatioDecider,
            deconv_class:type[InfoGANGenerator], layer_decider:LayerDecider,
            wiremask_bbo:bool, norm_wiremask:bool, input_partner_die:bool, input_alignment_mask:bool, 
            input_next_block:int, use_ready_layers_mask:bool, use_alignment_constraint:bool, set_vision_to_zero:int, set_canvas_to_zero:int,
            input_adjacent_terminal_mask:bool, input_adjacent_block_mask:bool, use_adjacent_terminal_constraint:bool, use_adjacent_block_constraint:bool,
            input_thermal_mask:bool, input_power_mask:bool,
        ) -> None:
        super().__init__()
        n_grid = int(num_act ** 0.5)
        self.shared_encoder = shared_encoder # global encoder for the whole chip
        self.transposed_conv = deconv_class(n_grid, self.shared_encoder.final_shape) # deconv the output of shared_encoder to the same shape as the original canvas
        self.local_encoder = local_encoder # relatively light-weighted encoder for local convolution
        self.fusion_module = nn.Conv2d(2, 1, 1) # this module should output final score mask
        self.actor_update_shared_encoder = actor_update_shared_encoder
    
        self.register_buffer("redundancy", torch.tensor(0.0))
        self.wiremask_bbo = wiremask_bbo
        self.ratio_decider = ratio_decider
        self.ratio_range = ratio_decider.ratio_range if ratio_decider is not None else None

        self.layer_decider = layer_decider
        self.layer_decider_forward = layer_decider is not None
        self.norm_wiremask = norm_wiremask
        self.input_partner_die = input_partner_die
        self.input_alignment_mask = input_alignment_mask
        self.input_next_block = input_next_block
        self.use_ready_layers_mask = use_ready_layers_mask
        self.use_alignment_constraint = use_alignment_constraint
        self.set_vision_to_zero = set_vision_to_zero
        self.set_canvas_to_zero = set_canvas_to_zero
        self.input_adjacent_terminal_mask = input_adjacent_terminal_mask
        self.input_adjacent_block_mask = input_adjacent_block_mask
        self.use_adjacent_terminal_constraint = use_adjacent_terminal_constraint
        self.use_adjacent_block_constraint = use_adjacent_block_constraint
        self.input_thermal_mask = input_thermal_mask
        self.input_power_mask = input_power_mask



    def get_device(self) -> torch.device:
        return self.redundancy.device

    
    def forward(self, obs:Batch, state:Union[dict, Batch, np.ndarray]=None, info=None) -> Tuple[Batch, torch.Tensor]:
        """
        Input:
            obs: the real observation, is a Batch.
            state: the hidden state of the model, is a torch.Tensor.
        Output:
            probs: the probability of each action, is a torch.Tensor.
            hidden: other info, such as other_probs dist.
        """
        device = self.get_device()
        canvas = obs["canvas"].to(device) # [B, D, H, W], D is the number of layers
        wiremask = obs["wiremask"].to(device) # [B, H, W]
        position_mask = obs["position_mask"].to(device) # [B, H, W]
        wiremask_next = obs["wiremask_next"].to(device) # [B, H, W]
        position_mask_next = obs["position_mask_next"].to(device) # [B, H, W]
        grid_area_next = torch.from_numpy(obs["grid_area_next"]).to(device=device, dtype=wiremask.dtype) # [B]

        if self.set_canvas_to_zero:
            canvas = torch.zeros_like(canvas)

        if self.norm_wiremask:
            num_grid = wiremask.shape[-1]
            num_net = torch.from_numpy(obs["num_net"]).to(device=device, dtype=wiremask.dtype) # [B]
            wiremask = wiremask / num_net
            wiremask_next = wiremask_next / num_net

        if not self.wiremask_bbo:
            if self.input_partner_die:
                num_die = canvas.shape[1]
                stacked_mask = [canvas[:,d] for d in range(num_die)] + [wiremask, position_mask] # [B, D+2, H, W]
                stacked_mask = torch.stack(stacked_mask, dim=1)
            else:
                layer_idx = obs["layer_idx"] # [B]
                b = np.arange(canvas.shape[0]) # [B]
                stacked_mask = torch.stack([canvas[b,layer_idx], wiremask, position_mask], dim=1) # [B, 3, H, W]
            
            if self.input_next_block == 1:
                stacked_mask = torch.cat([stacked_mask, wiremask_next.unsqueeze(1), position_mask_next.unsqueeze(1)], dim=1) # [B, C+2, H, W]
            
            if self.input_alignment_mask:
                alignment_mask = obs["alignment_mask"].to(device) # [B, H, W]
                stacked_mask = torch.cat([stacked_mask, alignment_mask.unsqueeze(1)], dim=1) # [B, C+1, H, W]

            if self.input_adjacent_terminal_mask:
                adjacent_terminal_mask = obs["adjacent_terminal_mask"].to(device) # [B, H, W]
                stacked_mask = torch.cat([stacked_mask, adjacent_terminal_mask.unsqueeze(1)], dim=1) # [B, C+3, H, W]

            if self.input_adjacent_block_mask:
                adjacent_block_mask = obs["adjacent_block_mask"].to(device)
                stacked_mask = torch.cat([stacked_mask, adjacent_block_mask.unsqueeze(1)], dim=1) # [B, C+4, H, W]

            if self.input_thermal_mask:
                thermal_mask = obs["thermal_mask"].to(device) # [B, D, H, W]
                stacked_mask = torch.cat([stacked_mask, thermal_mask], dim=1)
            
            if self.input_power_mask:
                power_mask = obs["power_mask"].to(device) # [B, D, H, W]
                stacked_mask = torch.cat([stacked_mask, power_mask], dim=1)

            assert stacked_mask.shape[1] == self.shared_encoder.input_channels, f"stacked_mask.shape[1] = {stacked_mask.shape[1]}, self.shared_encoder.input_channels = {self.shared_encoder.input_channels}"

            # set vision to zero, ablation study
            if self.set_vision_to_zero:
                stacked_mask = torch.zeros_like(stacked_mask)

            with torch.set_grad_enabled(self.actor_update_shared_encoder):
                global_feat = self.shared_encoder(stacked_mask, to_torch_as(obs['graph_data'], stacked_mask) if self.shared_encoder.graph else None) # [B, C]

            if self.ratio_decider is not None:
                ratio_mean, ratio_std = self.ratio_decider(global_feat, stacked_mask, grid_area_next) # [B]

            if self.layer_decider is not None and self.layer_decider_forward:
                layer_idx = torch.from_numpy(obs["layer_idx"]).to(device=device, dtype=torch.long) # [B]
                num_blk_without_placing_order = torch.from_numpy(obs["num_blk_without_placing_order"]).to(device=device, dtype=stacked_mask.dtype) # [B, num_die]
                
                # sequence feature 
                if self.layer_decider.async_place_input_sequence:
                    seq_feat = obs["sequence_feature"]
                else:
                    seq_feat = None

                # layer sequence
                if self.layer_decider.input_layer_sequence:
                    layer_sequence = torch.from_numpy(obs["layer_sequence"]).to(device=device, dtype=stacked_mask.dtype) # [B, ep_len]
                    layer_sequence_mask = torch.from_numpy(obs["layer_sequence_mask"]).to(device=device, dtype=torch.bool) # [B, ep_len]
                    layer_sequence_len = torch.from_numpy(obs["layer_sequence_len"]).to(device=device, dtype=torch.long) # [B]
                else:
                    layer_sequence = layer_sequence_mask = layer_sequence_len = None

                layer_prob = self.layer_decider(global_feat, stacked_mask, layer_idx, num_blk_without_placing_order, seq_feat, layer_sequence, layer_sequence_mask, layer_sequence_len) # [B, num_layer]
                
                if self.use_ready_layers_mask:
                    ready_layers = torch.from_numpy(obs["ready_layers"]).to(device=device, dtype=layer_prob.dtype) # [B, num_layer]
                    layer_prob = torch.where(ready_layers > 0, layer_prob, -1e8)
                layer_prob = torch.softmax(layer_prob, dim=-1)
                layer_decider_forward_finished = True
            else:
                layer_decider_forward_finished = False


                
            deconv_out = self.transposed_conv(global_feat) # [b,c,h,w]
            local_encoder_out = self.local_encoder(stacked_mask) # [b,c,h,w]
            merge_input = torch.cat([deconv_out, local_encoder_out], dim=1)
            score = self.fusion_module(merge_input)
            score = rearrange(score, 'b 1 h w -> b (h w)')

        else:
            layer_decider_forward_finished = False
            score = -1 * rearrange(wiremask, 'b h w -> b (h w)')
        

        # available_mask to decide final available position
        binary_alignment_mask = obs["binary_alignment_mask"].to(device) if self.use_alignment_constraint else torch.zeros_like(position_mask) # [B, H, W]
        position_mask_loose = obs["position_mask_loose"].to(device) # [B, H, W]
        adjacent_terminal_mask = obs["adjacent_terminal_mask"].to(device) if self.use_adjacent_terminal_constraint else torch.zeros_like(position_mask) # [B, H, W]
        binary_adjacent_block_mask = obs["binary_adjacent_block_mask"].to(device) if self.use_adjacent_block_constraint else torch.zeros_like(position_mask) # [B, H, W]
        boundary_mask = obs["boundary_mask"].to(device) # [B, H, W]

        # only 0+0 == 0
        available_mask = position_mask + binary_alignment_mask + adjacent_terminal_mask + binary_adjacent_block_mask

        # drop adjacent block mask constraint
        num_avail_pos = (available_mask == 0).sum(dim=(1,2))
        no_avail_pos = torch.where(num_avail_pos == 0)[0]
        if len(no_avail_pos) > 0:
            available_mask[no_avail_pos] = (position_mask + binary_alignment_mask + adjacent_terminal_mask)[no_avail_pos]

        # drop adjacent terminal mask constraint
        num_avail_pos = (available_mask == 0).sum(dim=(1,2)) # [B]
        no_avail_pos = torch.where(num_avail_pos == 0)[0]
        if len(no_avail_pos) > 0:
            available_mask[no_avail_pos] = (position_mask + binary_alignment_mask)[no_avail_pos]
        
        # drop alignment mask constraint
        num_avail_pos = (available_mask == 0).sum(dim=(1,2)) # [B]
        no_avail_pos = torch.where(num_avail_pos == 0)[0]
        if len(no_avail_pos) > 0:
            available_mask[no_avail_pos] = position_mask[no_avail_pos]

        # drop tight position mask constraint
        num_avail_pos = (available_mask == 0).sum(dim=(1,2)) # [B]
        no_avail_pos = torch.where(num_avail_pos == 0)[0]
        if len(no_avail_pos) > 0:
            available_mask[no_avail_pos] = position_mask_loose[no_avail_pos]

        # drop all constraints, except boundary mask
        num_avail_pos = (available_mask == 0).sum(dim=(1,2)) # [B]
        no_avail_pos = torch.where(num_avail_pos == 0)[0]
        if len(no_avail_pos) > 0:
            available_mask[no_avail_pos] = boundary_mask[no_avail_pos] # only position_mask is still not satisfied, use boundary_mask


        scoremask = score + rearrange(available_mask, 'b h w -> b (h w)') * -1e12
        probs_pos = torch.softmax(scoremask, dim=-1)
        hidden = None

        probs = Batch(pos=probs_pos)
        if self.ratio_decider is not None:
            probs["ratio"] = Batch(mean=ratio_mean, std=ratio_std)
        if layer_decider_forward_finished:
            probs["layer"] = layer_prob
            
        return probs, hidden
