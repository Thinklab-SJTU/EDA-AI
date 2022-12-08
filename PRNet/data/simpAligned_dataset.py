import os
from data.base_dataset import BaseDataset, get_params, get_transform
from data.image_folder import make_dataset
from PIL import Image
import torch
from tqdm import trange


class SimpAlignedDataset(BaseDataset):
    """A dataset class for paired image dataset.

    It assumes that the directory '/path/to/data/train' contains image pairs in the form of {A,B}.
    During test time, you need to prepare a directory '/path/to/data/test'.
    """

    def __init__(self, opt):
        """Initialize this dataset class.

        Parameters:
            opt (Option class) -- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        BaseDataset.__init__(self, opt)
        # input_folder = ['adaptec2_input', 'adaptec5_input', 'bigblue1_input', 'bigblue2_input', 'newblue1_input']
        # real_folder = ['adaptec2', 'adaptec5', 'bigblue1', 'bigblue2', 'newblue1']
        input_folder = ['bigblue1_obstacleInput_more_64']
        real_folder = ['bigblue1_obstacle_more_64']
        # input_folder = ['bigblue1_obstacleInput_more_64', 'newblue1_obstacleInput_more_64']
        # real_folder = ['bigblue1_obstacle_more_64', 'newblue1_obstacle_more_64']
        # input_folder = ['bigblue1_obstacleInput_64']
        # real_folder = ['bigblue1_obstacle_64']

        self.dir_A = [os.path.join(opt.dataroot, folder) for folder in input_folder]  # get the image directory
        self.dir_B = [os.path.join(opt.dataroot, folder) for folder in real_folder]
        self.A_paths = sorted(make_dataset(self.dir_A, opt.max_dataset_size))[:8000]  # get image paths
        self.B_paths = sorted(make_dataset(self.dir_B, opt.max_dataset_size))[:8000]
        self.input_nc = self.opt.input_nc
        self.output_nc = self.opt.output_nc

        # first execution
        # self.list_A = []  # list
        # self.list_B = []  # list
        #
        # for idx in trange(len(self.A_paths)):
        #     A_path = self.A_paths[idx]
        #     B_path = self.B_paths[idx]
        #     A = Image.open(A_path)
        #     B = Image.open(B_path).convert('L')
        #     # A = A.rotate(270)
        #     # B = B.rotate(270)
        #
        #     transform_params = get_params(self.opt, B)
        #     A_transform = get_transform(self.opt, transform_params, grayscale=(self.input_nc == 1))
        #     B_transform = get_transform(self.opt, transform_params, grayscale=(self.output_nc == 1))
        #
        #     A = A_transform(A)
        #     B = B_transform(B)
        #
        #     self.list_A.append(A)
        #     self.list_B.append(B)

        # torch.save({"list_A": self.list_A, "list_B": self.list_B}, "ispd08/obstacle_64_multi.pt")
        # torch.save({"list_A": self.list_A, "list_B": self.list_B}, "ispd08/obstacle_more_64.pt")
        # torch.save({"list_A": self.list_A, "list_B": self.list_B}, "ispd08/obstacle_64.pt")
        # torch.save({"list_A": self.list_A, "list_B": self.list_B}, "ispd08/obstacle_64_simple.pt")

        # further execution
        # ckpt = torch.load("ispd08/obstacle_more_64.pt")
        # ckpt = torch.load("ispd08/obstacle_64.pt")
        ckpt = torch.load("ispd08/obstacle_64_simple.pt")
        self.list_A = ckpt["list_A"]
        self.list_B = ckpt["list_B"]

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

        return {'A': self.list_A[index], 'B': self.list_B[index],
                'A_paths': self.A_paths[index], 'B_paths': self.B_paths[index]}

    def __len__(self):
        """Return the total number of images in the dataset."""
        return len(self.A_paths)
