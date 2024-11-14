import torch
from torch import nn
import torchvision
import numpy as np
from utils.utils import is_power_of_2
from .transformer import Transformer
from tianshou.data import Batch

class SharedEncoder(nn.Module):
    """
    This is the shared encoder, implemented by ResNet18.
    """
    def __init__(self, input_channels:int, graph:int, final_shape:tuple[int,int,int]):
        """final_shape: [C,H,W], the output can be reshaped to this shape."""
        super().__init__()
        assert final_shape[1] == final_shape[2], f"final_shape[1] != final_shape[2], got {final_shape[1]} and {final_shape[2]}"
        assert is_power_of_2(final_shape[1]) and is_power_of_2(final_shape[2]), f"final_shape[1] = {final_shape[1]}, final_shape[2] = {final_shape[2]}"
        self.encoder = torchvision.models.resnet18(weights='IMAGENET1K_V1')

        # change the first layer to accept input_channels
        self.input_channels = input_channels
        self.encoder.conv1 = nn.Conv2d(input_channels, 64, kernel_size=(7, 7), stride=(2, 2), padding=(3, 3), bias=False)

        # change the last layer to output final_shape
        self.final_shape = final_shape
        self.last_channel = np.prod(final_shape)
        self.encoder.fc = nn.Linear(512, self.last_channel)

        # graph
        self.graph = graph
        if graph:
            in_dim = 7
            hid_dim = 32
            if graph == 1:
                out_dim = self.last_channel
            elif graph == 2:
                out_dim = 32
                self.fc_graph = nn.Sequential(
                    nn.Linear(out_dim*2+self.last_channel, self.last_channel),
                    nn.ReLU(),
                    nn.Linear(self.last_channel, self.last_channel),
                )
            elif graph == 3:
                out_dim = 32
                self.fc_graph = nn.Sequential(
                    nn.Linear(out_dim+self.last_channel, self.last_channel),
                    nn.ReLU(),
                    nn.Linear(self.last_channel, self.last_channel),
                )
            else:
                raise ValueError(f"graph should be 0, 1, 2, 3, got {graph}")
            self.transformer = Transformer(in_dim, hid_dim, out_dim, n_layers=2)


    def forward(self, stacked_mask: torch.Tensor, graph_data_batch:Batch) -> torch.Tensor:
        """
        Input:
            stacked_mask: the stacked mask, is a torch.Tensor, shape = [B, C, H, W].
        Output:
            output: the output of the encoder, is a torch.Tensor, shape = [B, C].
        """
        input = stacked_mask
        output = self.encoder(input)

        if self.graph:
            graph_data = []
            for k in graph_data_batch.keys():
                if k in ['x', 'y', 'z', 'w', 'h', 'area', 'placed']:
                    graph_data.append(graph_data_batch[k])
            graph_data = torch.stack(graph_data, dim=-1).to(input.device)

            adj_mat = graph_data_batch['adj_mat_mov'].to(input.device)
            curr_node = graph_data_batch['idx'].long()
            order = graph_data_batch['order'].long() if 'order' in graph_data_batch.keys() else None

            emb_global, emb_local = self.transformer(graph_data, adj_mat, curr_node, order)
            if self.graph == 1:
                output = output + emb_global + emb_local
            elif self.graph == 2:
                output = torch.cat([output, emb_global, emb_local], dim=-1)
                output = self.fc_graph(output)
            elif self.graph == 3:
                output = torch.cat([output, emb_local+emb_global], dim=-1)
                output = self.fc_graph(output)
            else:
                raise ValueError(f"graph should be 0, 1, 2, 3, got {self.graph}")


        return output
    
# test
if __name__ == "__main__":
    input_channels = 3
    final_shape = (16, 8, 8)
    batch_size = 4
    encoder = SharedEncoder(input_channels, final_shape)
    x = torch.rand((batch_size, input_channels, 224, 224))
    print(encoder(x).shape)
