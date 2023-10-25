from torch.utils.data import Dataset
import os
from PIL import Image
from functools import reduce
import torch
from utils import get_transform, make_dataset, route2vertex

class AlignedDataset(Dataset):
    """A dataset class for paired image dataset.
    
    Ensure that images of routes and conditions both exist in 'HubRouter/dataset/', for example:
    
    HubRouter/
    |-- dataset/
    |   |-- adaptec1_64/
    |   |-- adaptec1_input_64/
    
    """

    def __init__(self, data_folder=None, condition_folder=None, route_folder=None, max_dataset_size=10000, restrict_c_channels=None, restrict_dir=None):
        """Initialize this dataset class.

        Parameters:
            data_folder: dataset/
            condition_folder: adaptec1_64/
            route_folder: adaptec1_input_64/
        """
        self.restrict_c_channels = restrict_c_channels
        if restrict_dir is not None:
            self.dir_condition = [os.path.join(data_folder, x, restrict_dir) for x in condition_folder.split(',')]
            self.dir_route = [os.path.join(data_folder, x, restrict_dir) for x in route_folder.split(',')]
        else:
            condition_folder = [os.path.join(data_folder, x) for x in condition_folder.split(',')]
            route_folder = [os.path.join(data_folder, x) for x in route_folder.split(',')]

            self.dir_condition = reduce(lambda x, y: x+y, [[os.path.join(dirs, dir) for dir in os.listdir(dirs)] for dirs in condition_folder])
            self.dir_route = reduce(lambda x, y: x+y, [[os.path.join(dirs, dir) for dir in os.listdir(dirs)] for dirs in route_folder])
        
        self.condition_paths = sorted(make_dataset(self.dir_condition, max_dataset_size))  # get image paths
        self.route_paths = sorted(make_dataset(self.dir_route, max_dataset_size))
        self.condition_transform = get_transform(grayscale=False)
        self.route_transform = get_transform(grayscale=True)

    def __getitem__(self, index):
        """Return a data point and its metadata information.

        Parameters:
            index - - a random integer for data indexing

        Returns a dictionary that contains condition, route, condition_paths and route_paths
            condition (tensor) - - an image in the condition domain
            route (tensor) - - its corresponding image in the target domain
            condition_paths (str) - - image paths
            route_paths (str) - - image paths (same as condition_paths)
        """

        condition_path = self.condition_paths[index]
        route_path = self.route_paths[index]
        
        assert condition_path.split('_')[-1] == route_path.split('_')[-1]
        condition = Image.open(condition_path).convert('RGB')
        route = Image.open(route_path).convert('L')

        condition = self.condition_transform(condition)
        if self.restrict_c_channels is not None:
            condition = condition[:1, :, :]
        route = self.route_transform(route)
        route_vertex = route2vertex(route)
        # route_stripe \in [0,1]
        route_stripe = (torch.sign((route_vertex+1).max(axis=1).values+(route_vertex+1).max(axis=2).values.T)).view(1,route_vertex.shape[1], route_vertex.shape[2])
        
        route = (route+1)/6 + (route_vertex+1)/6 + route_stripe*4/3 - 1

        return {'condition': condition, 'route': route, 'condition_paths': condition_path, 'route_paths': route_path}

    def __len__(self):
        """Return the total number of images in the dataset."""
        return len(self.condition_paths)