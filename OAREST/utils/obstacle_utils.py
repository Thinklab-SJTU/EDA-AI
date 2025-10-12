import numpy as np
import torch

def horizon_line_in_obstacles(line, y, input_obstacles):
    overlap = 0
    for i in range(len(input_obstacles)):
        if input_obstacles[i][1] < y < input_obstacles[i][3] and \
            (not (line[0] >= input_obstacles[i][2] or line[1] <= input_obstacles[i][0])):
            overlap += 1
    return overlap

def verticle_line_in_obstacles(line, x, input_obstacles):
    overlap = 0
    for i in range(len(input_obstacles)):
        if input_obstacles[i][0] < x < input_obstacles[i][2] and \
            (not (line[0] >= input_obstacles[i][3] or line[1] <= input_obstacles[i][1])):
            overlap += 1
    return overlap

def rec_overlap(existed_rec, rec):
    # Check if any points are inside the rectangle
    x_min, y_min, x_max, y_max = rec
    for i in range(len(existed_rec)):
        if min(existed_rec[i][2], x_max) >= max(existed_rec[i][0], x_min) and \
            min(existed_rec[i][3], y_max) >= max(existed_rec[i][1], y_min):
            return True
    return False

def points_in_rec(points, rec):
    # Check if the rectangle overlaps with any existing obstacles
    x_min, y_min, x_max, y_max = rec
    for i in range(len(points)):
        if x_min < points[i][0] < x_max and y_min < points[i][1] < y_max:
            return True
    return False

def count_overlap(input, input_obstacles, output):
    # input = new_input
    # input_obstacles = input_obstacle_batch[i]
    # output = output_list
    x_low, y_low, x_high, y_high = [], [], [], []
    overlap = 0
    for i in range(len(input)):
        x_low.append(input[i][0])
        y_low.append(input[i][1])
        x_high.append(input[i][0])
        y_high.append(input[i][1])
    for i in range(int(len(output) / 2)):
        x_idx = output[2*i]
        y_idx = output[2*i+1]
        x = input[x_idx][0]
        y = input[y_idx][1]
        y_low[x_idx] = min(y_low[x_idx], y)
        y_high[x_idx] = max(y_high[x_idx], y)
        x_low[y_idx] = min(x_low[y_idx], x)
        x_high[y_idx] = max(x_high[y_idx], x)
    for i in range(len(x_low)):
        if x_low[i] != x_high[i]:
            line = [x_low[i], x_high[i]]
            overlap += horizon_line_in_obstacles(line, input[i][1], input_obstacles)
            # continue # only count once for a rectilinear edge
        if y_low[i] != y_high[i]:
            line = [y_low[i], y_high[i]]
            overlap += verticle_line_in_obstacles(line, input[i][0], input_obstacles)
    return overlap

# def rand_tangle():
#     rec = np.random.rand(4)
#     rec[0], rec[2] = min(rec[0], rec[2]), max(rec[0], rec[2])
#     rec[1], rec[3] = min(rec[1], rec[3]), max(rec[1], rec[3])
#     return rec

# def rand_obstacles(input_batch, min_obstacle, max_obstacle):
#     input_obstacles_batch = []
#     complex_flag = 0
#     for i in range(len(input_batch)):
#         obstacles = []
#         num_obstacle = np.random.randint(min_obstacle, max_obstacle+1)
#         for j in range(num_obstacle):
#             rec = rand_tangle(complex_flag)
#             k = 0
#             while points_in_rec(input_batch[i], rec) or rec_overlap(obstacles, rec):
#                 # points are not in the rectangle
#                 rec = rand_tangle(complex_flag)
#                 k += 1
#                 if k >= 5:
#                     complex_flag = 1
#             obstacles.append(rec)
#         for j in range(max_obstacle-num_obstacle):
#             obstacles.append(-np.ones(4))
#         input_obstacles_batch.append(obstacles)
#     input_obstacles_batch = np.array(input_obstacles_batch)
#     return input_obstacles_batch

def rand_tangle(complex_flag):
    """
    Generate a rectangle with the bottom-left corner at a random position and ensure
    the top-right corner stays within the [0, 1] range.
    """
    x1, y1 = np.random.rand(2)  # Random coordinates for the bottom-left corner
    # Adjust the maximum width and height to ensure the top-right corner stays in [0, 1]
    x1, y1 = min(x1, 0.95), min(y1, 0.95)
    max_width = max(0.05, 1 - x1)  # Ensure at least 0.05 for width
    max_height = max(0.05, 1 - y1)  # Ensure at least 0.05 for height
    if complex_flag:
        width = 0.05
        height = 0.05
    else:
        width = np.random.uniform(0.05, max_width)  # Width range limited by available space
        height = np.random.uniform(0.05, max_height)  # Height range limited by available space
    x2, y2 = x1 + width, y1 + height  # Coordinates for the top-right corner
    return np.array([x1, y1, x2, y2])

def rand_obstacles(input_batch, min_obstacle, max_obstacle):
    """
    Generate a batch of random obstacle matrices.

    Args:
    - input_batch: Input points, each sample contains multiple points, shape (batch_size, num_points, 2).
    - min_obstacle: Minimum number of obstacles for each sample.
    - max_obstacle: Maximum number of obstacles for each sample.

    Returns:
    - input_obstacles_batch: A batch of randomly generated obstacles, shape (batch_size, max_obstacle, 4).
    """

    complex_flag = 0
    input_obstacles_batch = []
    for i in range(len(input_batch)):
        obstacles = []
        num_obstacle = np.random.randint(min_obstacle, max_obstacle + 1)  # Random number of obstacles
        for j in range(num_obstacle):
            rec = rand_tangle(complex_flag)
            k = 0
            # Retry rectangle generation to ensure points are outside the rectangle 
            # and there is no overlap with existing obstacles
            while points_in_rec(input_batch[i], rec) or rec_overlap(obstacles, rec):
                rec = rand_tangle(complex_flag)
                k += 1
                if k >= 5:
                    complex_flag = 1  # Trigger complex mode (placeholder for future logic)
            obstacles.append(rec)
        
        # Fill the remaining slots with placeholder `-np.ones(4)`
        for j in range(max_obstacle - num_obstacle):
            obstacles.append(-np.ones(4))
        
        input_obstacles_batch.append(obstacles)
    
    # Convert to NumPy array to standardize the batch output format
    input_obstacles_batch = np.array(input_obstacles_batch)
    return input_obstacles_batch

def obstacle_to_points(input_obstacle_batch):
    input_tensor = torch.tensor(input_obstacle_batch, dtype=torch.float32)  # (batch, obstacles, 4)
    points1 = input_tensor[..., :2]  # (batch, obstacles, 2)
    points2 = input_tensor[..., 2:]  # (batch, obstacles, 2)
    obstacle_points = torch.cat([points1, points2], dim=1)  # (batch, obstacles*2, 2)
    return obstacle_points

def obstacle_to_all_points(input_obstacle_batch):
    input_tensor = torch.tensor(input_obstacle_batch, dtype=torch.float32)  # (batch, obstacles, 4)
    left_down = input_tensor[..., :2]  # (batch, obstacles, 2)
    left_up = input_tensor[..., [0,3]]  # (batch, obstacles, 2)
    right_up = input_tensor[..., 2:]  # (batch, obstacles, 2)
    right_down = input_tensor[..., [2,1]]  # (batch, obstacles, 2)
    stacked_points = torch.stack([left_down.unsqueeze(2), left_up.unsqueeze(2), right_up.unsqueeze(2), right_down.unsqueeze(2)], dim=2)  # (batch, obstacles, 4, 2)
    reshaped_points = stacked_points.view(input_tensor.shape[0], input_tensor.shape[1] * 4, 2)
    return reshaped_points

def special_flip(a, degree):
    tmp_pointindex_to_newindex_dict = dict()
    tmp_newindex_to_pointindex_dict = dict()
    index = degree
    for i in range(len(a)-1, -1, -1):
        if a[i] >= degree and tmp_pointindex_to_newindex_dict.get(a[i], 0) == 0:
            tmp_pointindex_to_newindex_dict[a[i]] = index
            tmp_newindex_to_pointindex_dict[index] = a[i]
            index += 1
        elif a[i] < degree:
            tmp_pointindex_to_newindex_dict[a[i]] = a[i]
            tmp_newindex_to_pointindex_dict[a[i]] = a[i]
    return np.flip([tmp_pointindex_to_newindex_dict[x] for x in a])