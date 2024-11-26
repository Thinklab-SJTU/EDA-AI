import torch
import torchvision
import numpy as np

from utils import is_power_of_2
from tianshou.data import Batch
from torch import nn

from . import graph_model as GraphModel


def create_shared_encoder_backbone(shared_encoder_cls, input_channels:int, num_grid:int, final_shape:tuple[int,int,int]):
    if shared_encoder_cls is SharedEncoderLight:
        encoder = SharedEncoderLight(input_channels, num_grid, final_shape)

    elif shared_encoder_cls is SharedEncoder:
        encoder = torchvision.models.resnet18(weights='IMAGENET1K_V1')
        # change the first layer to accept input_channels
        encoder.conv1 = nn.Conv2d(input_channels, 64, kernel_size=(7, 7), stride=(2, 2), padding=(3, 3), bias=False)
        # change the last layer to output final_shape
        encoder.fc = nn.Linear(512, np.prod(final_shape))

    else:
        raise NotImplementedError(f"shared_encoder_cls = {shared_encoder_cls} is not implemented")
    
    return encoder


class SharedEncoderLight(nn.Module):
    """
    This is the shared encoder, implemented by a simple CNN, instead of ResNet18.
    """
    def __init__(self, input_channels:int, num_grid:int, final_shape:tuple[int,int,int]):
        super().__init__()

        final_down_sample_size = 16
        down_scale = round(np.log2(num_grid / final_down_sample_size))
        num_layer = down_scale

        base_dim = 32
        all_in_channels = [base_dim * 2 ** i for i in range(num_layer)]
        all_output_channels = [base_dim * 2 ** (i+1) for i in range(num_layer)]

        # initial conv
        conv = [
            nn.Conv2d(input_channels, all_in_channels[0], kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(all_in_channels[0]),
            nn.GELU(),
        ]

        # intermediate conv
        for c_in, c_out in zip(all_in_channels, all_output_channels):
            conv.extend([
                nn.Conv2d(c_in, c_out, kernel_size=3, stride=2, padding=1),
                nn.BatchNorm2d(c_out),
                nn.GELU(),
                nn.MaxPool2d(kernel_size=2, stride=2),
            ])
        
        # final conv
        # conv.extend([
        #     nn.Conv2d(all_output_channels[-1], all_output_channels[-1], kernel_size=3, stride=1, padding=1),
        #     nn.AdaptiveAvgPool2d((1, 1)),
        # ])
        conv.extend([
            nn.Conv2d(all_output_channels[-1], all_output_channels[-1], kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(all_output_channels[-1]),
            nn.GELU(),
            nn.AdaptiveMaxPool2d((1, 1)),
        ])


        self.conv = nn.Sequential(*conv)

        # final fc
        output_channels = np.prod(final_shape)
        self.fc = nn.Sequential(
            nn.Linear(all_output_channels[-1], output_channels),
            nn.GELU(),
            nn.Linear(output_channels, output_channels),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.conv(x)
        x = x.view(x.size(0), -1)
        x = self.fc(x)
        return x



class SharedEncoder(nn.Module):
    """
    This is the shared encoder, implemented by ResNet18.
    """
    def __init__(self, input_channels:int, graph:int, final_shape:tuple[int,int,int], num_grid:int, shared_encoder_cls:type[SharedEncoderLight], episode_len:int):
        """
        final_shape: [C,H,W], the output can be reshaped to this shape.
        """
        super().__init__()
        assert final_shape[1] == final_shape[2], f"final_shape[1] != final_shape[2], got {final_shape[1]} and {final_shape[2]}"
        assert is_power_of_2(final_shape[1]) and is_power_of_2(final_shape[2]), f"final_shape[1] = {final_shape[1]}, final_shape[2] = {final_shape[2]}"
        assert is_power_of_2(num_grid), f"num_grid = {num_grid} is not power of 2"
        
        self.input_channels = input_channels
        self.final_shape = final_shape
        self.num_grid = num_grid
        self.last_channel = np.prod(final_shape)

        self.encoder = create_shared_encoder_backbone(shared_encoder_cls, input_channels, num_grid, final_shape)

        # graph
        self.graph = graph
        if graph > 0:
            self.transformer, self.fc_graph = GraphModel.create_graph_model(graph, self.last_channel, episode_len)

    def forward(self, stacked_mask: torch.Tensor, graph_data_batch:Batch) -> torch.Tensor:
        """
        Input:
            stacked_mask: the stacked mask, is a torch.Tensor, shape = [B, C, H, W].
        Output:
            output: the output of the encoder, is a torch.Tensor, shape = [B, C].
        """
        output = self.encoder(stacked_mask)
        if self.graph > 0:
            output = GraphModel.forward_graph_model(self.graph, self.transformer, self.fc_graph, output, graph_data_batch)

        return output