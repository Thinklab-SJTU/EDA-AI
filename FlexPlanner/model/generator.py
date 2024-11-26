import torch
from torch import nn
import numpy as np
from utils import is_power_of_2


class BaseGenerator(nn.Module):
    def __init__(self, n_grid:int, input_shape:tuple[int, int, int]) -> None:
        super().__init__()

        assert is_power_of_2(n_grid), f"n_grid = {n_grid}"
        assert is_power_of_2(input_shape[1]) and is_power_of_2(input_shape[2]), f"input_shape = {input_shape}"
        assert input_shape[1] == input_shape[2], f"input_shape[1] != input_shape[2], got {input_shape[1]} and {input_shape[2]}"

        self.n_grid = n_grid
        self.input_shape = input_shape
    

    def forward(self, x:torch.Tensor) -> torch.Tensor:
        """Input is from SharedEncoder, shape is [B,C]"""
        raise NotImplementedError
    

    def count_parameters(self) -> int:
        num_params = sum(p.numel() for p in self.parameters() if p.requires_grad)
        print(f"[INFO][{self.__class__.__name__}] number of parameters = {num_params}")
        return num_params


class TransposedConv(BaseGenerator):
    def __init__(self, n_grid:int, input_shape:tuple[int, int, int]) -> None:
        super().__init__(n_grid, input_shape)

        deconv = []
        upscale = round(np.log2(n_grid // input_shape[-1]))
        in_channels = list(reversed([2**(i+1) for i in range(upscale)]))
        in_channels[0] = input_shape[0]
        out_channels = list(reversed([2**i for i in range(upscale)]))


        for in_c, out_c in zip(in_channels, out_channels):
            deconv.append(nn.ConvTranspose2d(in_c, out_c, kernel_size=3, stride=2, padding=1, output_padding=1, dilation=1))
            deconv.append(nn.ReLU())
        deconv.pop(-1)
        self.deconv = nn.Sequential(*deconv)
        # self.count_parameters()


    def forward(self, x:torch.Tensor) -> torch.Tensor:
        """Input is from SharedEncoder, shape is [B,C]"""
        x = x.reshape((-1,) + self.input_shape)
        assert x.shape[1:] == self.input_shape, f"In TransposedConv, x.shape is {x.shape} but should be {(-1,) + self.input_shape}"
        x = self.deconv(x)
        return x
    

class InfoGANGenerator(BaseGenerator):
    def __init__(self, n_grid:int, input_shape:tuple[int, int, int]) -> None:
        super().__init__(n_grid, input_shape)
        input_dim = np.prod(input_shape)
        hidden_dim1 = 32
        hidden_dim2 = 16
        self.hidden_dim1 = hidden_dim1
        self.hidden_dim2 = hidden_dim2

        self.init_size = input_shape[-1]  # Initial size before upsampling
        self.l1 = nn.Sequential(nn.Linear(input_dim, hidden_dim1 * self.init_size ** 2))

        # init blocks
        conv_blocks = [
            nn.BatchNorm2d(hidden_dim1), 
            nn.Upsample(scale_factor=2)
        ]

        # intermediate blocks
        num_layer = int(np.log2(n_grid // self.init_size)) - 1
        for _ in range(num_layer):
            conv_blocks.extend([
                nn.Conv2d(hidden_dim1, hidden_dim1, 3, stride=1, padding=1),
                nn.BatchNorm2d(hidden_dim1, 0.8),
                nn.LeakyReLU(0.2, inplace=True),
                nn.Upsample(scale_factor=2),
            ])

        # final blocks
        conv_blocks.extend([
            nn.Conv2d(hidden_dim1, hidden_dim2, 3, stride=1, padding=1),
            nn.BatchNorm2d(hidden_dim2, 0.8),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Conv2d(hidden_dim2, 1, 3, stride=1, padding=1),
            nn.Tanh(),
        ])
        
        self.conv_blocks = nn.Sequential(*conv_blocks)

        # self.count_parameters()
    

    def forward(self, x:torch.Tensor) -> torch.Tensor:
        """Input is from SharedEncoder, shape is [B,C]"""
        gen_input = x
        out = self.l1(gen_input)
        out = out.view(out.shape[0], self.hidden_dim1, self.init_size, self.init_size)
        img = self.conv_blocks(out)
        return img
