from torch.utils.data import Dataset
import numpy as np
import os
import torch
from utils import rand_obstacles, obstacle_to_points, obstacle_to_all_points
from utils import read_data

class RSMTDataset(Dataset):

    def __init__(self, min_degree, max_degree, dimension, min_obstacle, max_obstacle):
        self.min_degree = min_degree
        self.max_degree = max_degree
        self.dimension = dimension
        self.min_obstacle = min_obstacle
        self.max_obstacle = max_obstacle

    def __getitem__(self, index):
        input_batch = np.random.rand(1, self.max_degree, self.dimension)
        num_pins = np.random.randint(self.min_degree, self.max_degree+1)
        if num_pins < self.max_degree:
            input_batch[0, num_pins:] = -1
        if self.max_obstacle == 0:
            input_obstacles_batch = - np.ones([input_batch.shape[0], 1, 4])
        else:
            input_obstacles_batch = rand_obstacles(input_batch, self.min_obstacle, self.max_obstacle)
        input_obstacle_points_batch = obstacle_to_all_points(input_obstacles_batch)

        return input_batch[0], input_obstacles_batch[0], input_obstacle_points_batch[0]

    def __len__(self):
        """Return the total number of samples in the dataset."""
        return 125000

class RESTDataset(Dataset):
    def __init__(self, degree, min_obstacle, max_obstacle):
        self.degree = degree
        self.min_obstacle = min_obstacle
        self.max_obstacle = max_obstacle
        self.input_batch = read_data(f"test_set/test{degree}.txt")

    def __getitem__(self, index):
        input_batch = self.input_batch[index:index+1]
        if self.max_obstacle == 0:
            input_obstacles_batch = - np.ones([input_batch.shape[0], 1, 4])
        else:
            input_obstacles_batch = rand_obstacles(input_batch, self.min_obstacle, self.max_obstacle)
        input_obstacle_points_batch = obstacle_to_all_points(input_obstacles_batch)

        return input_batch[0], input_obstacles_batch[0], input_obstacle_points_batch[0]

    def __len__(self):
        """Return the total number of samples in the dataset."""
        return len(self.input_batch)

class ISPDDataset(Dataset):
    def __init__(self, min_degree, max_degree):
        self.min_degree = min_degree
        self.max_degree = max_degree
        if os.path.exists(f"dataset/ispd_test_{min_degree}_{max_degree}.np"):
            self.input_batch = np.load(f"dataset/ispd_test_{min_degree}_{max_degree}.np")
        else:
            with open(f"dataset/ispd_test_{min_degree}_{max_degree}.int", 'r') as f:
                input_batch = f.readlines()
                # max: 4301400, min: -2000
                input_batch = [np.array(eval(x)) for x in input_batch] # ensure the input is positive
                input_batch = [(x - np.min(x)) / (np.max(x) - np.min(x)) if np.max(x) != np.min(x) else x for x in input_batch]
            
            # padding
            input_batch = np.array([
                np.pad(x, ((0, max_degree - len(x)), (0, 0)), constant_values=0.) 
                for x in input_batch
            ])
            np.save(f"dataset/ispd_test_{min_degree}_{max_degree}.np", input_batch)
            # [batch_size, max_degree, 2]
            self.input_batch = input_batch

    def __getitem__(self, index):
        input_batch = self.input_batch[index:index+1]
        input_obstacles_batch = - np.ones([input_batch.shape[0], 1, 4])
        input_obstacle_points_batch = obstacle_to_points(input_obstacles_batch)

        return input_batch[0], input_obstacles_batch[0], input_obstacle_points_batch[0]

    def __len__(self):
        """Return the total number of samples in the dataset."""
        return len(self.input_batch)