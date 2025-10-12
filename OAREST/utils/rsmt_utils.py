import numpy as np
import torch

def transform_inputs(inputs: torch.Tensor, t: int) -> torch.Tensor:
    """
    Transform input coordinates based on transformation index t.
    
    Args:
        inputs: Tensor of shape (batch_size, num_points, 2) containing x,y coordinates
        t: Transform index (0-7)
        
    Returns:
        Transformed coordinates tensor of same shape
    """
    # Split into x and y coordinates
    xs = inputs[..., 0]  # Using ... to maintain any leading dimensions
    ys = inputs[..., 1]

    # Swap x and y for t >= 4
    if t >= 4:
        xs, ys = ys, xs
    
    # Mirror x coordinates
    if t % 2 == 1:
        xs = 1 - xs
    # Mirror y coordinates
    if t % 4 >= 2:
        ys = 1 - ys
    
    # Handle special case where coordinates are close to 2.0
    xs = torch.where(xs == 2.0, -1.0, xs)
    ys = torch.where(ys == 2.0, -1.0, ys)
    
    # Stack back into original format
    return torch.stack([xs, ys], dim=-1)

def transform_obstacles(obstacles: torch.Tensor, t: int) -> torch.Tensor:
    """
    Transform obstacle coordinates based on transformation index t.
    
    Args:
        obstacles: Tensor of shape (batch_size, num_obstacles, 4) containing x1,y1,x2,y2 coordinates
        t: Transform index (0-7)
        
    Returns:
        Transformed obstacles tensor of same shape
    """
    # Split into individual coordinates
    xs1 = obstacles[..., 0]
    ys1 = obstacles[..., 1]
    xs2 = obstacles[..., 2]
    ys2 = obstacles[..., 3]
    
    # Swap x and y coordinates for t >= 4
    if t >= 4:
        xs1, ys1 = ys1, xs1
        xs2, ys2 = ys2, xs2
        
    # Mirror x coordinates
    if t % 2 == 1:
        xs1, xs2 = 1 - xs2, 1 - xs1
        
    # Mirror y coordinates
    if t % 4 >= 2:
        ys1, ys2 = 1 - ys2, 1 - ys1
        
    # Handle special case where coordinates are close to 2.0
    for coords in [xs1, xs2, ys1, ys2]:
        coords.masked_fill_(coords == 2.0, -1)
    
    # Stack back into original format
    return torch.stack([xs1, ys1, xs2, ys2], dim=-1)

# def transform_inputs(inputs, t):
#     # 0 <= t <= 7
#     xs = inputs[:,:,0]
#     ys = inputs[:,:,1]
#     if t >= 4:
#         temp = xs
#         xs = ys
#         ys = temp
#     if t % 2 == 1:
#         xs = 1 - xs
#     if t % 4 >= 2:
#         ys = 1 - ys
#     xs[np.isclose(xs, 2.0)] = -1
#     ys[np.isclose(ys, 2.0)] = -1
#     return np.stack([xs, ys], -1)
    
# def transform_obstacles(obstacles, t):
#     # 0 <= t <= 7
#     xs1 = obstacles[:,:,0]
#     ys1 = obstacles[:,:,1]
#     xs2 = obstacles[:,:,2]
#     ys2 = obstacles[:,:,3]
#     if t >= 4:
#         xs1, ys1 = ys1, xs1
#         xs2, ys2 = ys2, xs2
#     if t % 2 == 1:
#         xs1, xs2 = 1 - xs2, 1 - xs1
#     if t % 4 >= 2:
#         ys1, ys2 = 1 - ys2, 1 - ys1
#     xs1[np.isclose(xs1, 2.0)] = -1
#     xs2[np.isclose(xs2, 2.0)] = -1
#     ys1[np.isclose(ys1, 2.0)] = -1
#     ys2[np.isclose(ys2, 2.0)] = -1
#     return np.stack([xs1, ys1, xs2, ys2], -1)

# def read_data(test_file):
#     with open(test_file, 'r') as f:
#         test_data = []
#         lines = f.readlines()
#         for line in lines:
#             line = line.strip().split(' ')
#             data = [float(coord) for coord in line]
#             data = np.array(data).reshape([-1, 2])
#             test_data.append(data)
#     return np.array(test_data)

def read_data(test_file):
    with open(test_file, 'r') as f:
        test_data = []
        lines = f.readlines()
        for line in lines:
            line = line.strip().split(' ')
            data = [float(coord) for coord in line]
            data = torch.tensor(data, dtype=torch.float32).view(-1, 2)
            test_data.append(data)
    return torch.stack(test_data)