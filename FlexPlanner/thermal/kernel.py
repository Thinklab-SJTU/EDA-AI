import torch
from torch import nn
import math

class MirrorPadding(nn.Module):
    def __init__(self):
        super().__init__()
    
    def forward(self, x:torch.Tensor) -> torch.Tensor:
        """x.shape should be [B, C, H, W]"""
        bs = x.shape[0]
        c = x.shape[1]

        flip_h = torch.flip(x, [2])
        flip_w = torch.flip(x, [3])
        flip_h_w = torch.flip(flip_h, [3])

        ans = torch.cat([
            flip_h_w, flip_h, flip_h_w,
            flip_w, x, flip_w,
            flip_h_w, flip_h, flip_h_w
        ]).reshape(bs, c, 3*x.shape[2], 3*x.shape[3])
        return ans
    

class ThermalConv2d(nn.Module):
    def __init__(self, num_grid_x:int, num_grid_y:int, thermal_mask:torch.Tensor=None) -> None:
        super().__init__()
        self.num_grid_x = num_grid_x
        self.num_grid_y = num_grid_y

        kernel = self.construct_kernel(num_grid_x+1, num_grid_y+1) if thermal_mask is None else thermal_mask
        self.register_buffer('kernel', kernel)

        self.pad = nn.ReflectionPad2d((num_grid_x//2, num_grid_x//2, num_grid_y//2, num_grid_y//2))
        # self.pad = MirrorPadding()

    
    def construct_kernel(self, num_grid_x:int, num_grid_y:int) -> torch.Tensor:
        # construct the kernel
        power_map = torch.zeros((num_grid_x, num_grid_y))
        x_start = math.floor((num_grid_x-1)/2)
        x_end = math.ceil((num_grid_x+1)/2)
        y_start = math.floor((num_grid_y-1)/2)
        y_end = math.ceil((num_grid_y+1)/2)
        power_map[x_start:x_end, y_start:y_end] = 1

        sum_power = power_map.sum().item()
        max_power = power_map.max().item()
        # print('sum_power = {}, max_power = {}'.format(sum_power, max_power))

        u = torch.arange(0, num_grid_x).float().reshape(-1, 1)
        x = torch.arange(0, num_grid_x).float().reshape(1, -1)
        A = 1 * math.pi * u * x / num_grid_x
        A = A.cos()

        v = torch.arange(0, num_grid_y).float().reshape(-1, 1)
        y = torch.arange(0, num_grid_y).float().reshape(1, -1)
        B = 1 * math.pi * v * y / num_grid_y
        B = B.cos()

        a_uv = A @ power_map @ B.T
        a_uv = a_uv / num_grid_x / num_grid_y

        u = torch.arange(0, num_grid_x).float().reshape(-1, 1)
        v = torch.arange(0, num_grid_y).float().reshape(1, -1)
        wu = 1 * math.pi * u / num_grid_x
        wv = 1 * math.pi * v / num_grid_y
        wu2_plus_wv2 = wu**2 + wv**2
        wu2_plus_wv2[0, 0] = 1

        x = torch.arange(0, num_grid_x).float().reshape(-1, 1)
        u = torch.arange(0, num_grid_x).float().reshape(1, -1)
        A = 1 * math.pi * u * x / num_grid_x
        A = A.cos()

        y = torch.arange(0, num_grid_y).float().reshape(-1, 1)
        v = torch.arange(0, num_grid_y).float().reshape(1, -1)
        B = 1 * math.pi * v * y / num_grid_y
        B = B.cos()

        kernel = A @ (a_uv / wu2_plus_wv2) @ B.T
        kernel = kernel.unsqueeze(0).unsqueeze(0)
        return kernel
    
    def forward(self, power_map:torch.Tensor) -> torch.Tensor:
        """
        power_map: (B, C, H, W) tensor, C is number of layer
        """
        power_map = self.pad(power_map)
        # print('padded power_map shape:', power_map.shape)
        return torch.nn.functional.conv2d(power_map, self.kernel)


class ThermalConvPB2d(nn.Module):
    def __init__(self, num_grid_x:int, num_grid_y:int, kernel_impulse:torch.Tensor, kernel_uniform:torch.Tensor, ambient_temperture:float) -> None:
        """
        num_grid_x and num_grid_y should be even number.
        kernel_impulse.shape should be [num_layer, num_layer, num_grid_x+1, num_grid_y+1].
        K_{oi} = thermal mask on layer o with point heat source on layer i.
        kernel_uniform.shape should be [num_layer, num_grid_x, num_grid_y].
        """
        super().__init__()
        assert num_grid_x % 2 == 0
        assert num_grid_y % 2 == 0
        assert kernel_impulse.shape[0] == kernel_impulse.shape[1]
        num_layer = kernel_impulse.shape[0]
        assert kernel_impulse.shape == (num_layer, num_layer, num_grid_x+1, num_grid_y+1)
        assert kernel_uniform.shape == (num_layer, num_grid_x, num_grid_y)

        self.num_grid_x = num_grid_x
        self.num_grid_y = num_grid_y
        self.num_layer = num_layer
        self.ambient_temperture = ambient_temperture

        kernel_impulse -= self.ambient_temperture
        kernel_uniform -= self.ambient_temperture

        self.register_buffer("kernel_impulse", kernel_impulse)
        self.register_buffer("kernel_uniform", kernel_uniform)

        self.pad = nn.ReflectionPad2d((num_grid_x//2, num_grid_x//2, num_grid_y//2, num_grid_y//2))

        intrinsic_error = self.construct_intrinsic_error()
        self.register_buffer("intrinsic_error", intrinsic_error)
        # print('intrinsic_error shape:', intrinsic_error.shape)
    
    @property
    def device(self) -> torch.device:
        return self.kernel_impulse.device
    

    def forward_dcmi(self, power_map:torch.Tensor) -> torch.Tensor:
        """
        DCMI = direct convolution with method of image.
        T0 = M00 * P0 + M01 * P1, where M_ij is the thermal mask obtained on layer i with point heat source on layer j.
        """
        power_map = self.pad(power_map)
        return torch.nn.functional.conv2d(power_map, self.kernel_impulse)
    

    def construct_intrinsic_error(self) -> torch.Tensor:
        power = 1.0 / self.num_grid_x / self.num_grid_y
        uniform_power_map = torch.full((1, self.num_layer, self.num_grid_x, self.num_grid_y), power).to(self.device)
        heatmap_dcmi = self.forward_dcmi(uniform_power_map)
        intrinsic_error = (heatmap_dcmi - self.kernel_uniform.unsqueeze(0)) / self.kernel_uniform.unsqueeze(0)
        return intrinsic_error

    def forward(self, power_map:torch.Tensor) -> torch.Tensor:
        """
        power_map: (B, C, H, W) tensor, C is number of layer
        """
        heatmap = self.forward_dcmi(power_map)
        heatmap = heatmap / (1 + self.intrinsic_error) + self.ambient_temperture
        return heatmap


if __name__ == '__main__':
    import os
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
    import matplotlib.pyplot as plt
    import numpy as np
    from thermal import thermal_fft_func
    plt.switch_backend('agg')

    num_grid_x = num_grid_y = 128
    num_layer = 2
    batch_size = 2
    ambient_temperture = 318.15

    kernel_impulse = torch.from_numpy(np.load(f'thermal/hotspot/impulse-grid={num_grid_x+1}.npy')).to(torch.float32)
    kernel_uniform = torch.from_numpy(np.load(f'thermal/hotspot/uniform-grid={num_grid_x}.npy')).to(torch.float32)

    conv = ThermalConvPB2d(num_grid_x, num_grid_y, kernel_impulse, kernel_uniform, ambient_temperture)

    # power_map = torch.rand(batch_size, num_layer, num_grid_x, num_grid_y)

    # power_map = torch.ones(batch_size, num_layer, num_grid_x, num_grid_y) / num_grid_x / num_grid_y

    power_map = torch.zeros(batch_size, num_layer, num_grid_x, num_grid_y)
    for r in [0.1, 0.3, 0.5, 0.7, 0.9]:
        w = 0.1 * num_grid_x
        xs = int(num_grid_x * r) - int(w/2)
        xe = xs + int(w)
        h = 0.1 * num_grid_y
        ys = int(num_grid_y * r) - int(h/2)
        ye = ys + int(h)
        power_map[:, :, xs:xe, ys:ye] = 1 / 10

    # power_map = torch.zeros(batch_size, num_layer, num_grid_x, num_grid_y)
    # xs = math.floor((num_grid_x-1)/2)
    # xe = math.ceil((num_grid_x+1)/2)
    # ys = math.floor((num_grid_y-1)/2)
    # ye = math.ceil((num_grid_y+1)/2)
    # power_map[:, 0, xs:xe, ys:ye] = 1 / (xe-xs) / (ye-ys)

    heatmap = conv(power_map)
    print('heatmap shape:', heatmap.shape)

    heatmap_dcmi = conv.forward_dcmi(power_map)

    heatmap_fft = thermal_fft_func(power_map) + conv.ambient_temperture

    nrow = 7
    fig, axes = plt.subplots(nrow, num_layer, figsize=(5*num_layer, 5*nrow))
    for i in range(num_layer):
        # impulse kernel
        r_idx = 0
        im = axes[r_idx, i].imshow(kernel_impulse[0, i].detach().numpy(), interpolation='none')
        temp_max = kernel_impulse[0, i].max().item()
        temp_min = kernel_impulse[0, i].min().item()
        axes[r_idx, i].set_title(f'impulse kernel layer {i}\nmax={temp_max:.2f}, min={temp_min:.2f}')
        fig.colorbar(im, ax=axes[r_idx, i])
        r_idx += 1

        # uniform kernel
        im = axes[r_idx, i].imshow(kernel_uniform[i].detach().numpy(), interpolation='none')
        temp_max = kernel_uniform[i].max().item()
        temp_min = kernel_uniform[i].min().item()
        axes[r_idx, i].set_title(f'uniform kernel layer {i}\nmax={temp_max:.2f}, min={temp_min:.2f}')
        fig.colorbar(im, ax=axes[r_idx, i])
        r_idx += 1

        # intrinsic error
        im = axes[r_idx, i].imshow(conv.intrinsic_error[0, i].detach().numpy(), interpolation='none')
        axes[r_idx, i].set_title(f'intrinsic error layer {i}')
        fig.colorbar(im, ax=axes[r_idx, i])
        r_idx += 1

        # power map
        im = axes[r_idx, i].imshow(power_map[0, i].detach().numpy(), interpolation='none')
        axes[r_idx, i].set_title(f'power map layer {i}')
        fig.colorbar(im, ax=axes[r_idx, i])
        r_idx += 1

        # heatmap_dcmi
        im = axes[r_idx, i].imshow(heatmap_dcmi[0, i].detach().numpy(), interpolation='none')
        temp_max = heatmap_dcmi[0, i].max().item()
        temp_min = heatmap_dcmi[0, i].min().item()
        axes[r_idx, i].set_title(f'heat map dcmi layer {i}\nmax={temp_max:.2f}, min={temp_min:.2f}')
        fig.colorbar(im, ax=axes[r_idx, i])
        r_idx += 1

        # heatmap
        im = axes[r_idx, i].imshow(heatmap[0, i].detach().numpy(), interpolation='none')
        temp_max = heatmap[0, i].max().item()
        temp_min = heatmap[0, i].min().item()
        axes[r_idx, i].set_title(f'heat map layer {i}\nmax={temp_max:.2f}, min={temp_min:.2f}')
        fig.colorbar(im, ax=axes[r_idx, i])
        r_idx += 1

        # heatmap_fft
        im = axes[r_idx, i].imshow(heatmap_fft[0, i].detach().numpy(), interpolation='none')
        temp_max = heatmap_fft[0, i].max().item()
        temp_min = heatmap_fft[0, i].min().item()
        axes[r_idx, i].set_title(f'heat map fft layer {i}\nmax={temp_max:.2f}, min={temp_min:.2f}')
        fig.colorbar(im, ax=axes[r_idx, i])


    plt.savefig('thermal/heatmap.png')
    plt.close()
