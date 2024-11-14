import torch
import math
from einops import rearrange

class ThermalFFT(torch.autograd.Function):
    @staticmethod
    def forward(ctx, power_map: torch.Tensor):
        """
        Use Poisson's equation and DCT to calculate temperature profile
        Input: power_map: (batch, channel, num_grid_x, num_grid_y)
        Output: temperature_profile: (batch, channel, num_grid_x, num_grid_y)
        """
        num_grid_x = power_map.shape[2]
        num_grid_y = power_map.shape[3]
        u = torch.arange(0, num_grid_x).float().reshape(-1, 1)
        x = torch.arange(0, num_grid_x).float().reshape(1, -1)
        A = 1 * math.pi * u * x / num_grid_x
        A = A.cos()

        v = torch.arange(0, num_grid_y).float().reshape(-1, 1)
        y = torch.arange(0, num_grid_y).float().reshape(1, -1)
        B = 1 * math.pi * v * y / num_grid_y
        B = B.cos()

        A = rearrange(A, 'u x -> 1 1 u x')
        B = rearrange(B, 'v y -> 1 1 y v') # transpose B

        # with shape (batch, channel, num_grid_x, num_grid_y)
        a_uv = A @ power_map @ B
        a_uv = a_uv / num_grid_x / num_grid_y

        u = torch.arange(0, num_grid_x).float().reshape(-1, 1)
        v = torch.arange(0, num_grid_y).float().reshape(1, -1)
        wu = 1 * math.pi * u / num_grid_x
        wv = 1 * math.pi * v / num_grid_y
        wu2_plus_wv2 = wu**2 + wv**2
        wu2_plus_wv2[0, 0] = 1
        wu2_plus_wv2 = rearrange(wu2_plus_wv2, 'u v -> 1 1 u v')

        x = torch.arange(0, num_grid_x).float().reshape(-1, 1)
        u = torch.arange(0, num_grid_x).float().reshape(1, -1)
        A = 1 * math.pi * u * x / num_grid_x
        A = A.cos()

        y = torch.arange(0, num_grid_y).float().reshape(-1, 1)
        v = torch.arange(0, num_grid_y).float().reshape(1, -1)
        B = 1 * math.pi * v * y / num_grid_y
        B = B.cos()

        A = rearrange(A, 'x u -> 1 1 x u')
        B = rearrange(B, 'y v -> 1 1 v y') # transpose B

        temperature_profile = A @ (a_uv / wu2_plus_wv2) @ B
        return temperature_profile


    @staticmethod
    def backward(ctx, grad_temperature_profile: torch.Tensor):
        return None

def thermal_fft_func(power_map: torch.Tensor):
    """
    Use Poisson's equation and DCT to calculate temperature profile
    Input: power_map: (batch, channel, num_grid_x, num_grid_y)
    Output: temperature_profile: (batch, channel, num_grid_x, num_grid_y)
    """
    return ThermalFFT.apply(power_map)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    plt.switch_backend('agg')

    num_grid_x = 128
    num_grid_y = 128

    bs = 3

    power_map = torch.zeros((bs, 1, num_grid_x, num_grid_y))
    for bs in range(power_map.shape[0]):
        # random set each sample
        if bs == 0:
            power_map[bs, 0, 20:30, 20:30] = 1
            power_map[bs, 0, 50:60, 50:60] = 1
            power_map[bs, 0, 80:90, 80:90] = 1
        elif bs == 1:
            power_map[bs, 0, 20:30, 20:30] = 1
            power_map[bs, 0, 50:60, 50:60] = 1
            power_map[bs, 0, 80:90, 80:90] = 1
            power_map[bs, 0, 110:120, 110:120] = 1
        elif bs == 2:
            power_map[bs, 0, 20:30, 20:30] = 1
            power_map[bs, 0, 50:60, 50:60] = 1
            power_map[bs, 0, 80:90, 80:90] = 1
            power_map[bs, 0, 110:120, 110:120] = 1
            power_map[bs, 0, 40:50, 40:50] = 1
            power_map[bs, 0, 70:80, 70:80] = 1
            power_map[bs, 0, 100:110, 100:110] = 1

    for bs in range(power_map.shape[0]):
        plt.imshow(power_map[bs,:].squeeze().numpy(), cmap='hot')
        plt.colorbar()
        plt.title('power map, bs = {}'.format(bs))
        plt.savefig(f'power_map_{bs}.png')
        plt.close()

    result = thermal_fft_func(power_map)
    for bs in range(result.shape[0]):
        plt.imshow(result[bs,0], cmap='hot')
        plt.colorbar()
        plt.title('temperature map, bs = {}'.format(bs))
        plt.savefig(f'result_{bs}.png')
        plt.close()


    # concat above imgs, use cv2
    power_map_imgs_paths = [f'power_map_{bs}.png' for bs in range(power_map.shape[0])]
    result_imgs_paths = [f'result_{bs}.png' for bs in range(result.shape[0])]
    import cv2
    power_map_imgs = [cv2.imread(img) for img in power_map_imgs_paths]
    result_imgs = [cv2.imread(img) for img in result_imgs_paths]
    power_map_imgs = np.concatenate(power_map_imgs, axis=1)
    result_imgs = np.concatenate(result_imgs, axis=1)
    # concat two imgs
    imgs = np.concatenate([power_map_imgs, result_imgs], axis=0)
    # save
    cv2.imwrite('thermal_fft.png', imgs)
    # remove
    import os
    for img in power_map_imgs_paths + result_imgs_paths:
        os.remove(img)



        
