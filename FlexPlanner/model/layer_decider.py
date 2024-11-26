import torch
from torch import nn
from einops.layers.torch import Rearrange
from tianshou.data import Batch, to_torch_as
from torch.nn.utils.rnn import pack_padded_sequence, pad_packed_sequence
from . import sequence_encoder as SeqEnc

class LayerDecider(nn.Module):
    """
    In asynchronous layer decision, the LayerDecider is responsible for deciding the layer to access the next block to place.
    """
    def __init__(self, num_die:int, input_dim:int, hidden_dim:int, share_with_critics:bool, 
                 die_embedding:bool, async_place_input_sequence:str, input_layer_sequence:bool):
        """
        share_with_critics: if True, use the input from shared encoder to decide the layer.
        async_place_input_sequence: if not None, input will be [batch_size, num_die, num_blk_each_die, channels].
        input_layer_sequence: if True, input will be [batch_size, num_total_blk, 1].
        """
        super().__init__()
        self.num_die = num_die
        self.enable_die_embedding = die_embedding
        self.async_place_input_sequence = async_place_input_sequence
        self.share_with_critics = share_with_critics
        self.input_layer_sequence = input_layer_sequence


        # CNN backbone
        if share_with_critics:
            print("[INFO] Use shared encoder in LayerDecider")
            self.layer_decider = nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim//2),
                nn.ReLU(),
                nn.Linear(hidden_dim//2, hidden_dim//4)
            )
            input_dim_final_layer = hidden_dim // 4
        else:
            self.layer_decider = nn.Sequential(
                nn.Conv2d(input_dim, hidden_dim, 3, 1, 1),
                nn.ReLU(),
                nn.MaxPool2d(2, 2),
                nn.Conv2d(hidden_dim, hidden_dim*2, 3, 1, 1),
                nn.ReLU(),
                nn.MaxPool2d(2, 2),
                nn.Conv2d(hidden_dim*2, hidden_dim*4, 3, 1, 1),
                nn.AdaptiveAvgPool2d((1, 1)),
                Rearrange('b c 1 1 -> b c'),
            )
            input_dim_final_layer = hidden_dim * 4

        
        # die embedding
        if self.enable_die_embedding:
            self.die_embedding = SeqEnc.DieEmbedding(num_die, num_die, trainable=False)
            input_dim_final_layer += 2*num_die # die_embedding, dim = num_die; num_blk_without_placing_order, dim = num_die.
            print("[INFO] Use die_embedding in LayerDecider")
            
        
        # multi-die place order sequence
        if self.async_place_input_sequence:
            seq_in_channels_each_die = 8
            seq_out_channels_each_die = 4
            input_dim_final_layer += num_die * seq_out_channels_each_die
            
            if async_place_input_sequence == 'CNN':
                self.seq_encoder = SeqEnc.SequenceEncoderCNN(seq_in_channels_each_die, 8, seq_out_channels_each_die, num_die)
            elif async_place_input_sequence == 'CNN2':
                self.seq_encoder = SeqEnc.SequenceEncoderCNN2(seq_in_channels_each_die, 8, seq_out_channels_each_die, num_die)
            elif async_place_input_sequence == 'Transformer2':
                self.seq_encoder = SeqEnc.Transformer2(seq_in_channels_each_die, 16, seq_out_channels_each_die, num_layers=2)
            else:
                raise NotImplementedError(f"async_place_input_sequence = {async_place_input_sequence} is not implemented")
            print(f"[INFO] Use sequence encoder {self.seq_encoder.__class__.__name__} in LayerDecider")
        

        # layer sequence
        if self.input_layer_sequence:
            assert self.enable_die_embedding, "input_layer_sequence requires die_embedding"
            rnn_in_size = num_die
            rnn_hid_size = num_die*2
            input_dim_final_layer += rnn_hid_size
            self.rnn = nn.GRU(rnn_in_size, rnn_hid_size, num_layers=1, batch_first=True)
            print("[INFO] Use layer sequence in LayerDecider")

        # final header
        if input_dim_final_layer == num_die:
            self.final_header = nn.Identity() # no additional feature to concat, so the final header is just an identity
        else:
            self.final_header = nn.Sequential(
                nn.Linear(input_dim_final_layer, input_dim_final_layer),
                nn.ReLU(),
                nn.Linear(input_dim_final_layer, num_die)
            )

    def rnn_forward(self, layer_sequence:torch.LongTensor, layer_sequence_mask:torch.BoolTensor, layer_sequence_len:torch.LongTensor) -> torch.Tensor:
        """return the final step of rnn output."""
        layer_sequence = layer_sequence.long() # [B, ep_len]
        layer_sequence = torch.where(layer_sequence_mask == 1, self.num_die, layer_sequence) # 1 means invalid
        layer_sequence = self.die_embedding(layer_sequence) # change int to embedding, [B, ep_len] -> [B, ep_len, C]
        layer_sequence = pack_padded_sequence(layer_sequence, layer_sequence_len.cpu(), batch_first=True, enforce_sorted=False)

        rnn_out, _ = self.rnn(layer_sequence)
        rnn_out, lens_unpacked = pad_packed_sequence(rnn_out, batch_first=True)
        rnn_out = rnn_out[torch.arange(rnn_out.size(0)), lens_unpacked-1]
        return rnn_out


    def forward(self, output_of_shared_encoder:torch.Tensor, stacked_mask:torch.Tensor, layer_idx:torch.LongTensor, 
                num_blk_without_placing_order:torch.Tensor, sequence_feature:Batch, 
                layer_sequence:torch.Tensor, layer_sequence_mask:torch.BoolTensor, layer_sequence_len:torch.LongTensor) -> torch.Tensor:
        """
        output_of_shared_encoder & stacked_mask: [B, C, H, W]. These features are vision modality.
        sequence_feature: [B, num_die, num_blk_each_die, C]. These features are sequence modality.
        """
        x = output_of_shared_encoder if self.share_with_critics else stacked_mask
        x = self.layer_decider(x) # [B, C]

        if self.enable_die_embedding:
            die_embedding = self.die_embedding(layer_idx) # [B, num_die]
            x = torch.cat([x, die_embedding, num_blk_without_placing_order], dim=-1)
        
        if self.async_place_input_sequence:
            # [B, num_die, num_blk, C] -> [B, num_die]
            sequence_feature = self.seq_encoder(to_torch_as(sequence_feature, x))
            x = torch.cat([x, sequence_feature], dim=-1)
        
        if self.input_layer_sequence:
            rnn_out_final_step = self.rnn_forward(layer_sequence, layer_sequence_mask, layer_sequence_len)
            x = torch.cat([x, rnn_out_final_step], dim=-1) # concat last step
            
        x = self.final_header(x)
        return x

