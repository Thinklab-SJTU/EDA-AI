from torch.utils.data import Dataset
import os
from PIL import Image
from functools import reduce
import torch
from utils import get_transform, make_dataset, route2vertex

class ConditionDataset(Dataset):
    """A dataset class for paired image dataset.
    
    Ensure that images of routes and conditions both exist in 'HubRouter/dataset/', for example:
    
    HubRouter/
    |-- dataset/
    |   |-- adaptec1_64/
    |   |-- adaptec1_input_64/
    
    """

    def __init__(self, data_folder=None, condition_folder=None, max_dataset_size=10000):
        """Initialize this dataset class.

        Parameters:
            data_folder: dataset/
            condition_folder: adaptec1_64/
        """
        
        condition_folder = [os.path.join(data_folder, x) for x in condition_folder.split(',')]
        
        self.dir_condition = reduce(lambda x, y: x+y, [[os.path.join(dirs, dir) for dir in os.listdir(dirs)] for dirs in condition_folder])
        
        self.condition_paths = sorted(make_dataset(self.dir_condition, max_dataset_size))  # get image paths
        self.condition_transform = get_transform(grayscale=False)

    def __getitem__(self, index):
        """Return a data point and its metadata information.

        Parameters:
            index - - a random integer for data indexing

        Returns a dictionary that contains condition, route, condition_paths and route_paths
            condition (tensor) - - an image in the condition domain
            condition_paths (str) - - image paths
        """

        condition_path = self.condition_paths[index]
        condition = Image.open(condition_path).convert('RGB')

        condition = self.condition_transform(condition)

        return {'condition': condition, 'condition_paths': condition_path}

    def __len__(self):
        """Return the total number of images in the dataset."""
        return len(self.condition_paths)