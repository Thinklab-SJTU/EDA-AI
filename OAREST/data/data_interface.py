import gc
import torch
from torch.utils.data import DataLoader, Dataset
import pytorch_lightning as pl

from utils import instantiate_from_config

class CachedDataset(Dataset):
    def __init__(self, cached_data):
        self.cached_data = cached_data

    def __getitem__(self, index):
        return self.cached_data[index]

    def __len__(self):
        return len(self.cached_data)

class DInterface(pl.LightningDataModule):
    def __init__(self, data_config, batch_size):
        super().__init__()
        self.num_workers = data_config.num_workers
        self.dataset_name = data_config.dataset
        self.batch_size = batch_size
        min_degree = data_config.train_val_data_config.params.min_degree
        max_degree = data_config.train_val_data_config.params.max_degree
        self.dataset = instantiate_from_config(data_config.train_val_data_config)
        self.split()
        if self.dataset_name == 'random':
            self.test_dataset = None
        elif self.dataset_name == 'rest' and min_degree == max_degree and max_degree in list(range(5, 51, 5)):
            self.test_dataset = instantiate_from_config(data_config.test_data_config)
        elif self.dataset_name == 'ispd':
            self.test_dataset = instantiate_from_config(data_config.ISPD_data_config)

    def setup(self, stage=None):
        # Assign train/val datasets for use in dataloaders
        self.trainset = self.train_dataloader()
        self.valset = self.val_dataloader()

        # Assign test dataset for use in dataloader(s)
        self.testset = self.test_dataloader()

    def train_dataloader(self):
        gc.collect()
        torch.cuda.empty_cache()
        return DataLoader(self.train_dataset,
                          batch_size=self.batch_size,
                          shuffle=True,
                          num_workers=self.num_workers, 
                          persistent_workers=True,
                          pin_memory=True,
                          prefetch_factor=2)
    
    def val_dataloader(self):
        # return DataLoader(self.val_dataset,
        #                   batch_size=self.batch_size,
        #                   shuffle=False,
        #                   num_workers=self.num_workers)
        val_dataset = CachedDataset(self.val_cached_data)
        return DataLoader(val_dataset, 
                          batch_size=self.batch_size, 
                          shuffle=False, 
                          num_workers=self.num_workers, 
                          persistent_workers=True,
                          pin_memory=True,
                          prefetch_factor=2)

    def test_dataloader(self):
        # return DataLoader(self.test_dataset,
        #                   batch_size=self.batch_size,
        #                   shuffle=False,
        #                   num_workers=self.num_workers)
        if self.test_dataset:
            return DataLoader(self.test_dataset,
                            batch_size=self.batch_size,
                            shuffle=False,
                            num_workers=self.num_workers, 
                            persistent_workers=True,
                            pin_memory=True,
                            prefetch_factor=2)  
        else:
            test_dataset = CachedDataset(self.test_cached_data)
            return DataLoader(test_dataset, 
                              batch_size=self.batch_size, 
                              shuffle=False, 
                              num_workers=self.num_workers, 
                              persistent_workers=True,
                              pin_memory=True,
                              prefetch_factor=2)    
                          
    def all_dataloader(self):
        return DataLoader(self.dataset,
                          batch_size=self.batch_size,
                          shuffle=False,
                          num_workers=self.num_workers, 
                          persistent_workers=True,
                          pin_memory=True,
                          prefetch_factor=2)

    def split(self, test=True):
        train_size = int(0.80 * len(self.dataset))
        val_size = int(0.10 * len(self.dataset))
        test_size = int(0.10 * len(self.dataset))
        self.train_dataset, val_dataset_temp, test_dataset_temp = torch.utils.data.random_split(self.dataset, 
                                [train_size, val_size, test_size], generator=torch.Generator().manual_seed(0))

        self.val_cached_data = [val_dataset_temp[i] for i in range(len(val_dataset_temp))]
        self.test_cached_data = [test_dataset_temp[i] for i in range(len(test_dataset_temp))]