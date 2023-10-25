from torch.utils.data import Dataset
import numpy as np

class RSMTDataset(Dataset):

    def __init__(self, degree, dimension):
        self.degree = degree
        self.dimension = dimension

    def __getitem__(self, index):
        input_batch = np.random.rand(self.degree, self.dimension)

        return input_batch

    def __len__(self):
        """Return the total number of images in the dataset."""
        return 256 * 10000