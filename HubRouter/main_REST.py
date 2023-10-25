# -*- coding: utf-8 -*-
from omegaconf import OmegaConf
import torch
import lightning.pytorch as pl
from lightning.pytorch import Trainer
from lightning.pytorch.loggers import TensorBoardLogger

from models import ACInterface
from data import DInterface
from utils import load_callbacks, examine_dir
import argparse

torch.set_float32_matmul_precision('high')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, default=None, help='seed')
    parser.add_argument('--degree', type=int, default=None, help='seed')
    return parser.parse_args()

def train(args):
    config = OmegaConf.load(f"configs/REST.yaml")
    seed = args.seed if args.seed is not None else config.base.seed
    pl.seed_everything(seed)
    
    batch_size = config.base.batch_size
    ckpt_dir = config.base.ckpt_dir
    result_dir = config.base.result_dir
    log_dir = config.base.log_dir
    ckpt_path = None if config.base.ckpt_path == 'None' else config.base.ckpt_path
    model_name = config.base.model_name
    degree = args.degree if args.degree is not None else config.base.degree
    dimension = config.base.dimension
    
    examine_dir(ckpt_dir)
    examine_dir(result_dir)
    examine_dir(log_dir)
    examine_dir(ckpt_dir+model_name+'_seed'+str(seed))

    data_module = DInterface(config.data, batch_size)
    
    model = ACInterface(config.model, 
                        model_name, 
                        result_dir, 
                        ckpt_dir, 
                        batch_size, 
                        degree, 
                        seed)
    callbacks = load_callbacks('val_lengths_epoch')

    logger = TensorBoardLogger(save_dir=log_dir, name=f'{model_name}_seed{seed}.log')
    trainer = Trainer(accelerator="gpu", 
                      max_epochs=config.base.max_epochs, 
                      callbacks=callbacks, 
                      logger=logger)
    if ckpt_path is not None:
        model.load_state_dict(torch.load(ckpt_path), strict=True)
    else:
        trainer.fit(model, data_module, ckpt_path=ckpt_path)
    
if __name__ == '__main__':
    args = get_args()
    train(args)

    
