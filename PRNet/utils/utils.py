"""This module contains simple helper functions """
from __future__ import print_function
import torch
import numpy as np
from PIL import Image
import os
from torch.autograd import Variable
import sys

np.random.seed(1)


def tensor2im(input_image, imtype=np.uint8):
    """"Converts a Tensor array into a numpy image array.

    Parameters:
        input_image (tensor) --  the input image tensor array
        imtype (type)        --  the desired type of the converted numpy array
    """
    if not isinstance(input_image, np.ndarray):
        if isinstance(input_image, torch.Tensor):  # get the data from a variable
            image_tensor = input_image.data
        else:
            return input_image
        image_numpy = image_tensor[0].cpu().float().numpy()  # convert it into a numpy array
        if image_numpy.shape[0] == 1:  # grayscale
            image_numpy = image_numpy[0]
            # image_numpy = np.heaviside(image_numpy, 0.) * 255.0
            image_numpy = (image_numpy + 1) / 2.0 * 255.0
            mode = 'L'
        else:
            image_numpy = (np.transpose(image_numpy,
                                        (1, 2, 0)) + 1) / 2.0 * 255.0  # post-processing: tranpose and scaling
            mode = 'RGB'
    else:  # if it is a numpy array, do nothing
        image_numpy = input_image
        mode = 'RGB'
    return image_numpy.astype(imtype), mode


def diagnose_network(net, name='network'):
    """Calculate and print the mean of average absolute(gradients)

    Parameters:
        net (torch network) -- Torch network
        name (str) -- the name of the network
    """
    mean = 0.0
    count = 0
    for param in net.parameters():
        if param.grad is not None:
            mean += torch.mean(torch.abs(param.grad.data))
            count += 1
    if count > 0:
        mean = mean / count
    print(name)
    print(mean)


def save_image(image_numpy, image_path, mode):
    """Save a numpy image to the disk

    Parameters:
        image_numpy (numpy array) -- input numpy array
        image_path (str)          -- the path of the image
    """

    image_pil = Image.fromarray(image_numpy, mode=mode)

    image_pil.save(image_path)


def print_numpy(x, val=True, shp=False):
    """Print the mean, min, max, median, std, and size of a numpy array

    Parameters:
        val (bool) -- if print the values of the numpy array
        shp (bool) -- if print the shape of the numpy array
    """
    x = x.astype(np.float64)
    if shp:
        print('shape,', x.shape)
    if val:
        x = x.flatten()
        print('mean = %3.3f, min = %3.3f, max = %3.3f, median = %3.3f, std=%3.3f' % (
            np.mean(x), np.min(x), np.max(x), np.median(x), np.std(x)))


def mkdirs(paths):
    """create empty directories if they don't exist

    Parameters:
        paths (str list) -- a list of directory paths
    """
    if isinstance(paths, list) and not isinstance(paths, str):
        for path in paths:
            mkdir(path)
    else:
        mkdir(paths)


def mkdir(path):
    """create a single empty directory if it didn't exist

    Parameters:
        path (str) -- a single directory path
    """
    if not os.path.exists(path):
        os.makedirs(path)


def make_folder(path, version):
    if not os.path.exists(os.path.join(path, version)):
        os.makedirs(os.path.join(path, version))


def tensor2var(x, grad=False):
    if torch.cuda.is_available():
        x = x.cuda()
    return Variable(x, requires_grad=grad)


def var2tensor(x):
    return x.data.cpu()


def var2numpy(x):
    return x.data.cpu().numpy()


def denorm(x):
    out = (x + 1) / 2
    return out.clamp_(0, 1)


def dfs(grid, r, c):
    grid[r][c] = 0
    n = len(grid)
    for x, y in [(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)]:
        if 0 <= x < n and 0 <= y < n and grid[x][y] > 0:
            dfs(grid, x, y)


def judge_constraint(grid_A, grid_B) -> bool:
    index = np.argwhere(grid_A > 0)

    if index.shape[-1] == 0:
        return True

    for j, k in zip(index[0], index[1]):
        if grid_B[j][k] < 0:
            return False

    # sys.setrecursionlimit(5000)
    try:
        dfs(grid_B, index[0][0], index[1][0])
    except RecursionError:
        return False

    for j, k in zip(index[0][1:], index[1][1:]):
        if grid_B[j][k] > 0:
            return False

    return True


def constraint_label(tensor_A, tensor_B):
    bs = tensor_A.shape[0]
    res = list()

    for i in range(bs):
        grid_A = tensor_A[i, 0, :, :].clone()
        grid_B = tensor_B[i, 0, :, :].clone()

        if judge_constraint(grid_A, grid_B):
            res.append(1.0)
        else:
            res.append(0.0)

    return torch.tensor(res)


def cal_jam(grid_A, grid_B):
    A_1 = grid_A[1]
    A_2 = grid_A[2]

    n = len(grid_B)

    zero_pad = torch.zeros([1, n])
    B_1 = torch.cat([grid_B[1:n], zero_pad], 0)
    B_1 = torch.gt(B_1, 0.) * torch.gt(grid_B, 0.)

    zero_pad = torch.zeros([n, 1])
    B_2 = torch.cat([grid_B[:, 1:n], zero_pad], 1)
    B_2 = torch.gt(B_2, 0.) * torch.gt(grid_B, 0.)

    res_1 = torch.ge((A_1 + 1.) / 2. * 255., B_1.float())
    res_2 = torch.ge((A_2 + 1.) / 2. * 255., B_2.float())

    return (np.argwhere(res_1 < 1).shape[1]) + (np.argwhere(res_2 < 1).shape[1])


def jam_label(tensor_A, tensor_B):
    bs = tensor_A.shape[0]
    res = list()

    for i in range(bs):
        grid_A = tensor_A[i, :, :, :].clone()
        grid_B = tensor_B[i, 0, :, :].clone()

        res.append(cal_jam(grid_A, grid_B))
    return torch.tensor(res).float()


def length_label(tensor_B):
    bs = tensor_B.shape[0]
    res = list()

    for i in range(bs):
        grid_B = tensor_B[i, 0, :, :]
        res.append(np.argwhere(grid_B > 0).shape[1])

    return torch.tensor(res).float()


def hpwl(tensor_A):
    bs = tensor_A.shape[0]
    res = list()

    for i in range(bs):
        pins = tensor_A[i, 0, :, :].clone()
        index = np.argwhere(pins > 0)
        arr_x, arr_y = index[0], index[1]
        min_x, max_x = torch.min(arr_x).item(), torch.max(arr_x).item()
        min_y, max_y = torch.min(arr_y).item(), torch.max(arr_y).item()
        x = max_x - min_x + 1
        y = max_y - min_y + 1
        res.append(x + y - 1)

    return torch.tensor(res).float()


