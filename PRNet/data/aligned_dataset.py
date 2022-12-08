import os
from data.base_dataset import BaseDataset, get_params, get_transform
from data.image_folder import make_dataset
from PIL import Image
import torch
from tqdm import trange


class AlignedDataset(BaseDataset):
    """A dataset class for paired image dataset.

    It assumes that the directories '/path/to/data/train/A' and '/path/to/data/train/B' contains images of A and B.
    During test time, you need to prepare directories '/path/to/data/test/A' and '/path/to/data/test/B'.
    """

    def __init__(self, opt):
        """Initialize this dataset class.

        Parameters:
            opt (Option class) -- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        BaseDataset.__init__(self, opt)

        input_folder = os.path.join(opt.dataroot, opt.A_folder)
        real_folder = os.path.join(opt.dataroot, opt.B_folder)

        self.dir_A = [os.path.join(opt.dataroot, folder) for folder in input_folder]  # get the image directory
        self.dir_B = [os.path.join(opt.dataroot, folder) for folder in real_folder]
        self.A_paths = sorted(make_dataset(self.dir_A, opt.max_dataset_size))  # get image paths
        self.B_paths = sorted(make_dataset(self.dir_B, opt.max_dataset_size))
        self.input_nc = self.opt.input_nc
        self.output_nc = self.opt.output_nc

    def __getitem__(self, index):
        """Return a data point and its metadata information.

        Parameters:
            index - - a random integer for data indexing

        Returns a dictionary that contains A, B, A_paths and B_paths
            A (tensor) - - an image in the input domain
            B (tensor) - - its corresponding image in the target domain
            A_paths (str) - - image paths
            B_paths (str) - - image paths (same as A_paths)
        """

        A_path = self.A_paths[idx]
        B_path = self.B_paths[idx]
        A = Image.open(A_path).convert('RGB')
        B = Image.open(B_path).convert('L')

        transform_params = get_params(self.opt, B)
        A_transform = get_transform(self.opt, transform_params, grayscale=(self.input_nc == 1))
        B_transform = get_transform(self.opt, transform_params, grayscale=(self.output_nc == 1))

        A = A_transform(A)
        B = B_transform(B)

        return {'A': A, 'B': B, 'A_paths': A_path, 'B_paths': B_path}

    def __len__(self):
        """Return the total number of images in the dataset."""
        return len(self.A_paths)
