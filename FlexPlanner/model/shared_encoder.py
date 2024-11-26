import torch
from torch import nn
import torchvision
import numpy as np
from utils.utils import is_power_of_2
from . import graph_model as GraphModel
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
        if graph > 0:
            self.transformer, self.fc_graph = GraphModel.create_graph_model(graph, self.last_channel)


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
            output = GraphModel.forward_graph_model(self.graph, self.transformer, self.fc_graph, output, graph_data_batch)

        return output
    
# test
if __name__ == "__main__":
    input_channels = 3
    final_shape = (16, 8, 8)
    batch_size = 4
    encoder = SharedEncoder(input_channels, final_shape)
    x = torch.rand((batch_size, input_channels, 224, 224))
    print(encoder(x).shape)
