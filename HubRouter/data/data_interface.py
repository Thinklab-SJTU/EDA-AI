import torch
from torch.utils.data import DataLoader
import lightning.pytorch as pl

from utils import instantiate_from_config
        
class DInterface(pl.LightningDataModule):
    def __init__(self, data_config, batch_size):
        super().__init__()
        self.num_workers = data_config.num_workers
        self.batch_size = batch_size
        self.dataset = instantiate_from_config(data_config)
        self.split()

    def setup(self, stage=None):
        # Assign train/val datasets for use in dataloaders
        self.trainset = self.train_dataloader()
        self.valset = self.val_dataloader()

        # Assign test dataset for use in dataloader(s)
        self.testset = self.test_dataloader()

    def train_dataloader(self):
        return DataLoader(self.train_dataset,
                          batch_size=self.batch_size,
                          shuffle=True,
                          num_workers=self.num_workers)

    def val_dataloader(self):
        return DataLoader(self.test_dataset,
                          batch_size=self.batch_size,
                          shuffle=False,
                          num_workers=self.num_workers)

    def test_dataloader(self):
        return DataLoader(self.test_dataset,
                          batch_size=self.batch_size,
                          shuffle=False,
                          num_workers=self.num_workers)
                          
    def all_dataloader(self):
        return DataLoader(self.dataset,
                          batch_size=self.batch_size,
                          shuffle=False,
                          num_workers=self.num_workers)             

    def split(self):
        train_size = int(0.8 * len(self.dataset))
        test_size = len(self.dataset) - train_size
        self.train_dataset, self.test_dataset = torch.utils.data.random_split(self.dataset, 
                                [train_size, test_size], generator=torch.Generator().manual_seed(0))