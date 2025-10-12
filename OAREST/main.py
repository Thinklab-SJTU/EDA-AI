# -*- coding: utf-8 -*-
import os
os.environ['CUDA_LAUNCH_BLOCKING'] = '1'

import glob
import torch
import argparse
import warnings
# Suppress all warnings
warnings.filterwarnings("ignore")
import pytorch_lightning as pl
from pytorch_lightning import Trainer
from pytorch_lightning.loggers import TensorBoardLogger
from omegaconf import OmegaConf

from models import ACInterface
from data import DInterface
from utils import load_callbacks, examine_dir

torch.set_float32_matmul_precision('high')
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', default="OAREST_test", help='config file')
    parser.add_argument('--seed', type=int, help='seed')
    parser.add_argument('--degree', type=int, help='degree')
    return parser.parse_args()

def train(args):
    config = OmegaConf.load(f"configs/{args.d}.yaml")
    seed = args.seed if args.seed is not None else config.base.seed
    pl.seed_everything(seed)
    
    mode = config.base.mode
    batch_size = config.base.batch_size
    img_dir = config.base.img_dir
    log_dir = config.base.log_dir
    ckpt_path = None if config.base.ckpt_path == 'None' else config.base.ckpt_path
    model_name = config.base.model_name
    min_degree = config.base.min_degree
    max_degree = args.degree if args.degree is not None else config.base.max_degree
    min_obstacle = config.base.min_obstacle
    max_obstacle = config.base.max_obstacle
    transformation = config.base.transformation
    assert max_degree >= min_degree >= 0
    assert max_obstacle >= min_obstacle >= 0
    # dimension = config.base.dimension
    obstacle_weight = config.base.obstacle_weight
    dataset = config.data.dataset

    config.model.actor_config.params.max_degree = max_degree
    config.model.critic_config.params.max_degree = max_degree
    config.data.train_val_data_config.params.min_degree = min_degree
    config.data.train_val_data_config.params.max_degree = max_degree
    config.data.train_val_data_config.params.min_obstacle = min_obstacle
    config.data.train_val_data_config.params.max_obstacle = max_obstacle
    if dataset == 'rest':
        config.data.test_data_config.params.degree = max_degree
        config.data.test_data_config.params.min_obstacle = min_obstacle
        config.data.test_data_config.params.max_obstacle = max_obstacle
    elif dataset == 'ispd':
        config.data.ISPD_data_config.params.max_degree = max_degree
        config.data.ISPD_data_config.params.min_degree = min_degree

    examine_dir(img_dir)
    examine_dir(log_dir)

    data_module = DInterface(config.data, batch_size)
    
    model = ACInterface(config.model, 
                        model_name, 
                        img_dir, 
                        batch_size, 
                        min_degree, 
                        max_degree, 
                        min_obstacle, 
                        max_obstacle, 
                        obstacle_weight, 
                        transformation, 
                        seed)
    callbacks = load_callbacks('val_metric')

    log_name = f'{model_name}_degree-[{min_degree},{max_degree}]_obstacle-[{min_obstacle},{max_obstacle}]_seed-{seed}.log'
    logger = TensorBoardLogger(save_dir=log_dir, name=log_name)
    root_log_dir = os.path.join(log_dir, log_name)

    trainer = Trainer(accelerator="gpu", 
                      min_epochs=config.base.max_epochs,
                      max_epochs=config.base.max_epochs, 
                      callbacks=callbacks, 
                      logger=logger)
    if ckpt_path is not None:
        if "DAC" in ckpt_path:
            # load states in REST
            checkpoint = torch.load(ckpt_path, map_location=model.device)
            # model.load_state_dict(checkpoint['state_dict'], strict=True)
            model.actor.load_state_dict(checkpoint['actor_state_dict'], strict=False)
            model.critic.load_state_dict(checkpoint['critic_state_dict'], strict=False)
        else:
            checkpoint = torch.load(ckpt_path, map_location=model.device)
            model.load_state_dict(checkpoint['state_dict'], strict=True)
    if mode in ['train', 'both']:
        trainer.fit(model, data_module)
    if mode in ['test', 'both']:
        trainer = Trainer(accelerator="gpu", callbacks=[], logger=False) # don't save checkpoints
        if ckpt_path is None:
            version_dirs = sorted(glob.glob(os.path.join(glob.escape(root_log_dir), "version_*")), key=os.path.getmtime)
            latest_version_dir = version_dirs[-1] if version_dirs else None
            checkpoint_dir = os.path.join(latest_version_dir, "checkpoints")
            best_checkpoint = glob.glob(f"{glob.escape(checkpoint_dir)}/best-epoch=*.ckpt")
            if best_checkpoint:
                best_checkpoint_path = best_checkpoint[0]
                print(f"Best checkpoint path: {best_checkpoint_path}")
                best_model = ACInterface.load_from_checkpoint(best_checkpoint_path)
                trainer.test(best_model, data_module)
            else:
                print("No best checkpoint found in the latest version.")
        else:
            trainer.test(model, data_module)

if __name__ == '__main__':
    args = get_args()
    train(args)
