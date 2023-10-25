# -*- coding: utf-8 -*-
from omegaconf import OmegaConf
import torch
import lightning.pytorch as pl

from models import DDPMInterface, GANInterface, VAEInterface, ACInterface
from data import DInterface
import os
import argparse
from time import time
import numpy as np
import pandas as pd
from utils import get_input_parameters_ibm, read_input_info
from REST_tool.REST_utils import Evaluator, transform_inputs
from preprocess.grid import Grid
import torchvision.transforms as transforms
import random
random.seed(333)

torch.set_float32_matmul_precision('high')

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', type=str, choices=['DDPM', 'VAE', 'GAN'], default='GAN', help='choose which model')
    parser.add_argument('--case', type=str, choices=['ISPD', 'REST', 'IBM', 'GRL8', 'GRL16'], default='IBM', help='case')
    parser.add_argument('--seed', type=int, default=None, help='seed')
    return parser.parse_args()

def inference(args):
    if args.case == 'IBM':
        ibm_flag = 1 
        args.case = 'ISPD'
    else:
        ibm_flag = 0
    if args.case[:3] == 'GRL':
        grl_flag = 1
        scale = args.case[3:]
        args.case = 'ISPD'
    else:
        grl_flag = 0
    config = OmegaConf.load(f"configs/FMCG-{args.model}-{args.case}.yaml")
    seed = args.seed if args.seed is not None else config.base.seed
    pl.seed_everything(seed)
    device = torch.device("cuda:0")
    
    batch_size = config.base.batch_size

    ckpt_dir = config.base.ckpt_dir
    result_dir = config.base.result_dir
    model_name = config.base.model_name
    if args.case == 'REST':
        data_config = OmegaConf.load(f"configs/test-{args.case}.yaml")
        data_module = DInterface(data_config.data, batch_size=batch_size)
        REST_model_config = OmegaConf.load("configs/REST.yaml")
        degree = data_config.data.params.restrict_dir.split('_')[-1]
        
        path = f'logs/REST_seed{seed}.log/'
        max_version = max([int(v.split('_')[-1]) for v in os.listdir(path)])
        # max_version = 0
        path = path + f'version_{max_version}/checkpoints/'
        rest_ckp_dir = path + list(filter(lambda x: x[:4] == 'best', os.listdir(path)))[0]
        # rest_ckp_dir = f'./REST_tool/pretrained/rsmt{degree}b.pt'
        REST_model = ACInterface(REST_model_config.model, 
                                 model_name, 
                                 result_dir, 
                                 ckpt_dir, 
                                 batch_size, 
                                 degree, 
                                 seed)
        checkpoint = torch.load(rest_ckp_dir, map_location=device)
        REST_model.actor.load_state_dict(checkpoint['actor_state_dict'])
        REST_model.actor.eval()
        
    evaluator = Evaluator(path=os.getcwd() + '/REST_tool/algorithms/libeval.so') # We borrow the evaluator from REST
    
    if args.model == 'DDPM':
        model = DDPMInterface(config.model, 
                              model_name, 
                              result_dir,
                              ckpt_dir, 
                              batch_size, 
                              seed)
    elif args.model == 'GAN':
        lambda_norm = config.base.lambda_norm
        lambda_focal = config.base.lambda_focal
        model = GANInterface(config.model, 
                             model_name, 
                             result_dir, 
                             ckpt_dir, 
                             batch_size, 
                             lambda_norm, 
                             lambda_focal, 
                             seed)
    elif args.model == 'VAE':
        latent_size = config.base.latent_size
        model = VAEInterface(config.model, 
                             model_name, 
                             result_dir, 
                             ckpt_dir, 
                             batch_size, 
                             seed, 
                             latent_size)
    
    path = f'logs/FMCG_{args.model}_{args.case}_64_seed{seed}.log/'
    max_version = max([int(v.split('_')[-1]) for v in os.listdir(path)])
    # max_version = 0
    path = path + f'version_{max_version}/checkpoints/'
    path = path + list(filter(lambda x: x[:4] == 'best', os.listdir(path)))[0]
    
    checkpoint = torch.load(path, map_location=device)
    model.load_state_dict(checkpoint["state_dict"])
    split_size = 64
    
    model = model.to(device)
    if args.model == 'DDPM':
        model.model.eval()
    else:
        model.eval()
    
    if ibm_flag == 1:
        benchmark = ['ibm01.modified.txt', 'ibm02.modified.txt', 'ibm03.modified.txt', 
             'ibm04.modified.txt', 'ibm05.modified.txt', 'ibm06.modified.txt']
        for i in range(len(benchmark)):
            scale = 64
            aligned_generated_routing_twl = 0.
            generation_time = 0.
            postprocess_time = 0.
            grid_info = read_input_info('preprocess/benchmark/' + benchmark[i])
            grid_param, nets = get_input_parameters_ibm(grid_info)
            env = Grid(grid_param)
            net_array = np.array([[net['net_name'], net['num_pin'], str(net['pins']), net['same_tile']] for net in nets.values()])
            net_pd = pd.DataFrame(net_array, columns = ['net_name', 'num_pin', 'pins', 'same_tile'])
            net_pd.num_pin = net_pd.num_pin.astype(int)
            net_pd.same_tile = net_pd.same_tile.apply(lambda x: eval(x))
            # net_pd_same_tile = net_pd[net_pd.same_tile == True]
            net_pd_not_same_tile = net_pd[net_pd.same_tile == False]
            net_pd_not_same_tile = net_pd_not_same_tile.sort_values('num_pin')

            env.planar_cap = torch.tensor(env.planar_cap).to(device)
            with torch.no_grad():
                for k in range(len(net_pd_not_same_tile)):
                    pins = eval(net_pd_not_same_tile.iloc[k, 2])
                    pins = np.array(pins)
                    condition_img = torch.zeros((1, 3, scale, scale)).to(device)
                    
                        
                    if max(pins[:, 0]) -  min(pins[:, 0]) <= 63 and max(pins[:, 1]) -  min(pins[:, 1]) <= 63:
                        random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0), min(pins[:, 0].min(), env.n_x-scale))
                        random_y = random.randint(max(pins[:, 1].max() - scale + 1, 0), min(pins[:, 1].min(), env.n_y-scale))
                        
                        condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x+scale, random_y:random_y+scale, 0]
                        condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x+scale, random_y:random_y+scale, 1]
                        for pin in pins:
                            condition_img[0, 0, pin[0]-random_x, pin[1]-random_y] = 255.
                        condition_img /= 255.
                        real_pin = condition_img[:, :1, :, :].detach()
                    else:
                        try:
                            random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0), min(pins[:, 0].min(), env.n_x-scale))
                        except:
                            random_x = min(pins[:, 0].min(), env.n_x-scale)
                        try:
                            random_y = random.randint(max(pins[:, 0].max() - scale + 1, 0), min(pins[:, 1].min(), env.n_y-scale))
                        except:
                            random_y = min(pins[:, 1].min(), env.n_y-scale)
                        
                        condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x+scale, random_y:random_y+scale, 0]
                        condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x+scale, random_y:random_y+scale, 1]
                        for pin in pins:
                            try:
                                condition_img[0, 0, pin[0]-random_x, pin[1]-random_y] = 255.
                            except:
                                pass
                        condition_img /= 255.
                        real_pin = condition_img[:, :1, :, :].detach()
                                            
                    t0 = time()
                    if args.model == 'DDPM':
                        # model.model.ddim_timesteps = 75
                        x_gen = model.model.ddim_sample(condition_img, guide_w=1.5)
                    elif args.model == 'GAN':
                        x_gen = model(condition_img)
                    elif args.model == 'VAE':
                        z = torch.randn(condition_img.shape[0], config.model.vae_config.params.latent_size).to(device)
                        x_gen = model.vae.inference(z, condition_img)
                    t1 = time()
        
                    generation_time += t1 - t0
                    x_mask = torch.sign((x_gen.sum(2).unsqueeze(2) + x_gen.sum(3).unsqueeze(-1)))
                    x_mask = (x_mask + 1)/2. # to [0,1]
                    x_gen = x_gen * x_mask
                    
                    # also add rectangular mask
                    for batch_j in range(real_pin.shape[0]):
                        pin_index = torch.where(real_pin[batch_j] == 1.)
                        mask_route = torch.zeros_like(real_pin[batch_j])
                        mask_route[0, pin_index[1].min() : pin_index[1].max() + 1, 
                               pin_index[2].min() : pin_index[2].max() + 1] =  1.
                        x_gen[batch_j] = x_gen[batch_j] * mask_route
                    
                    x_vertex = torch.max(x_gen, real_pin)
                    
                    t1 = time()
                    for test_case in x_vertex[:, 0, :, :]:
                        tc = test_case.detach().cpu().numpy()
                        if args.model in ['DDPM', 'VAE']:
                            tc = np.array(np.where(tc >= 5/6)).T
                        elif args.model in ['GAN']:
                            tc = np.array(np.where(tc >= 1/2)).T
                            
                        tc[:, 0] += random_x
                        tc[:, 1] += random_y
                        tc = np.concatenate((tc, pins))
                        tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                        
                        gst_length, sp_list, edge_list = evaluator.gst_rsmt(tc)
                        aligned_generated_routing_twl += gst_length

                        t2 = time()
                        postprocess_time += t2 - t1
                        
                        sp_list = np.int32(sp_list)
                        if len(sp_list) != 0:
                            points = np.concatenate([tc, sp_list], 0)
                        else:
                            points = tc

                        
                        for edge in range(len(edge_list)):
                            u = edge_list[edge][0]
                            v = edge_list[edge][1]
                            if points[u][0] == points[v][0]:
                                env.planar_cap[points[u][0], min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                            elif points[u][1] == points[v][1]:
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]), points[v][1], 0] -= 1
                            else:
                                env.planar_cap[points[u][0], min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]), points[v][1], 0] -= 1

                        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
                                f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                                f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                                f'total time: {generation_time + postprocess_time}', end='\r')
                print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
                f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                f'total time: {generation_time + postprocess_time}')
                
    elif grl_flag == 1:
        if scale == '8':
            benchmark = ['test_benchmark_88_' + str(i) + '.gr' for i in range(1, 11)]
        else:
            benchmark = ['test_benchmark_' + str(i) + '.gr' for i in range(1, 11)]
        
        scale = int(scale)
        total_overflow = 0.
        total_aligned_generated_routing_twl = 0.
        total_generation_time = 0.
        total_postprocess_time = 0.

        for i in range(len(benchmark)):
            aligned_generated_routing_twl = 0.
            generation_time = 0.
            postprocess_time = 0.
            grid_info = read_input_info('preprocess/benchmark/' + benchmark[i])
            grid_param, nets = get_input_parameters_ibm(grid_info)
            for i in range(len(grid_param['reduced_capacity'])):
                grid_param['reduced_capacity'][i][2] -= 1
                grid_param['reduced_capacity'][i][5] -= 1
            env = Grid(grid_param)
            net_array = np.array([[net['net_name'], net['num_pin'], str(net['pins']), net['same_tile']] for net in nets.values()])
            net_pd = pd.DataFrame(net_array, columns = ['net_name', 'num_pin', 'pins', 'same_tile'])
            net_pd.num_pin = net_pd.num_pin.astype(int)
            net_pd.same_tile = net_pd.same_tile.apply(lambda x: eval(x))
            # net_pd_same_tile = net_pd[net_pd.same_tile == True]
            net_pd_not_same_tile = net_pd[net_pd.same_tile == False]
            net_pd_not_same_tile = net_pd_not_same_tile.sort_values('num_pin')

            env.planar_cap = torch.tensor(env.planar_cap).to(device)
            with torch.no_grad():
                for k in range(len(net_pd_not_same_tile)):
                    
                    pins = eval(net_pd_not_same_tile.iloc[k, 2])
                    pins = np.array(pins)
                    condition_img = torch.zeros((1, 3, 64, 64)).to(device)
                    
                        
                    condition_img[0, 1, :scale, :scale] = env.planar_cap[:scale, :scale, 0]
                    condition_img[0, 2, :scale, :scale] = env.planar_cap[:scale, :scale, 1]
                    for pin in pins:
                        condition_img[0, 0, pin[0], pin[1]] = 255.
                    condition_img /= 255.
                    real_pin = condition_img[:, :1, :, :].detach()
                                            
                    t0 = time()
                    if args.model == 'DDPM':
                        # model.model.ddim_timesteps = 75
                        x_gen = model.model.ddim_sample(condition_img, guide_w=1.5)
                    elif args.model == 'GAN':
                        x_gen = model(condition_img)
                    elif args.model == 'VAE':
                        z = torch.randn(condition_img.shape[0], config.model.vae_config.params.latent_size).to(device)
                        x_gen = model.vae.inference(z, condition_img)
                    t1 = time()
        
                    generation_time += t1 - t0
                    x_mask = torch.sign((x_gen.sum(2).unsqueeze(2) + x_gen.sum(3).unsqueeze(-1)))
                    x_mask = (x_mask + 1)/2. # to [0,1]
                    x_gen = x_gen * x_mask
                    
                    # also add rectangular mask
                    for batch_j in range(real_pin.shape[0]):
                        pin_index = torch.where(real_pin[batch_j] == 1.)
                        mask_route = torch.zeros_like(real_pin[batch_j])
                        mask_route[0, pin_index[1].min() : pin_index[1].max() + 1, 
                               pin_index[2].min() : pin_index[2].max() + 1] =  1.
                        x_gen[batch_j] = x_gen[batch_j] * mask_route
                    
                    x_vertex = torch.max(x_gen, real_pin)
                    
                    for test_case in x_vertex[:, 0, :, :]:
                        t1 = time()
                        tc = test_case.detach().cpu().numpy()
                        if args.model in ['DDPM', 'VAE']:
                            tc = np.array(np.where(tc >= 5/6)).T
                        elif args.model in ['GAN']:
                            tc = np.array(np.where(tc >= 1/2)).T
                            
                        tc = np.concatenate((tc, pins))
                        tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                        
                        gst_length, sp_list, edge_list = evaluator.gst_rsmt(tc)
                        aligned_generated_routing_twl += gst_length

                        t2 = time()
                        postprocess_time += t2 - t1
                        
                        sp_list = np.int32(sp_list)
                        if len(sp_list) != 0:
                            points = np.concatenate([tc, sp_list], 0)
                        else:
                            points = tc

                        
                        for edge in range(len(edge_list)):
                            u = edge_list[edge][0]
                            v = edge_list[edge][1]
                            if points[u][0] == points[v][0]:
                                env.planar_cap[points[u][0], min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                            elif points[u][1] == points[v][1]:
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]), points[v][1], 0] -= 1
                            else:
                                env.planar_cap[points[u][0], min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]), points[v][1], 0] -= 1

                        print(f'benchmark: {benchmark[i]}, net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
                                f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                                f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                                f'total time: {generation_time + postprocess_time}', end='\r')
                print()
                total_overflow += (-env.planar_cap * (env.planar_cap < 0)).sum().item()
                total_aligned_generated_routing_twl += aligned_generated_routing_twl
                total_generation_time += generation_time
                total_postprocess_time += postprocess_time
        print(f'net: {k}/{len(net_pd_not_same_tile)}, total generated wl: {total_aligned_generated_routing_twl}, ' + \
                f'total generation time: {total_generation_time}s, total postprocess_time: {total_postprocess_time}s, ' + \
                f'total overflow: {total_overflow} ' + \
                f'total time: {total_generation_time + total_postprocess_time}')
            
    elif args.case == 'REST':
        aligned_real_routing_twl = 0.
        aligned_generated_routing_twl = 0.
        generation_time = 0.
        postprocess_time = 0.
        batch_id = 0
        num_batch = len(data_module.all_dataloader())
        print('degree:', degree)
        with torch.no_grad():
            for batch in data_module.all_dataloader():
                c = batch['condition'].to(device)
                x = batch['route'].to(device)
                
                real_pin = c[:, :1, :, :].detach()
                real_routing_length = torch.sum(x > 0.6)
                
                aligned_real_routing_twl += (real_routing_length.sum().item()-batch_size) / split_size
                
                t0 = time()
                if args.model == 'DDPM':
                    # model.model.ddim_timesteps = 75
                    x_gen = model.model.ddim_sample(c, guide_w=1.5)
                elif args.model == 'GAN':
                    x_gen = model(c)
                elif args.model == 'VAE':
                    z = torch.randn(c.shape[0], config.model.vae_config.params.latent_size).to(device)
                    x_gen = model.vae.inference(z, c)
                t1 = time()
    
                generation_time += t1 - t0
                x_mask = torch.sign((x_gen.sum(2).unsqueeze(2) + x_gen.sum(3).unsqueeze(-1))).to(device)
                x_mask = (x_mask + 1)/2. # to [0,1]
                x_gen = x_gen * x_mask
                
                x_vertex = torch.max(x_gen, real_pin)
                
                ## Geosteiner
                if data_config.data.rsmt == 'GEO':
                    for test_case in x_vertex[:, 0, :, :]:
                        tc = test_case.detach().cpu().numpy()
                        if args.model in ['DDPM', 'VAE']:
                            tc = np.array(np.where(tc >= 5/6)).T
                        elif args.model in ['GAN']:
                            tc = np.array(np.where(tc >= 1/2)).T
                        
                        gst_length, sp_list, edge_list = evaluator.gst_rsmt(tc)
                        aligned_generated_routing_twl += gst_length / split_size
                    t2 = time()
                    postprocess_time += t2 - t1
                
                ## REST
                elif data_config.data.rsmt == 'REST':                    
                    for b in range(len(x_vertex)):
                        tc = x_vertex[b, 0, :, :]
                        if args.model in ['DDPM']:
                            tc = torch.cat(torch.where(tc >= 5/6)).view(2, -1).T
                        elif args.model in ['VAE', 'GAN']:
                            tc = torch.cat(torch.where(tc >= 1/2)).view(2, -1).T
                        test_batch = tc.unsqueeze(0)
                        test_batch = test_batch / split_size
                        REST_model.actor.degree = test_batch.shape[1]
                        test_batch = test_batch.cpu().detach().numpy()
                        transformed_batches = []
                        for t in range(8):
                            transformed_batch = transform_inputs(test_batch, t)
                            transformed_batches.append(transformed_batch)
                        
                        transformed_batches = torch.tensor(np.array(transformed_batches)).squeeze(1).to(device)
                        
                        t1 = time()
                        outputs, _ =  REST_model.actor(transformed_batches, True)
                        t2 = time()
                        postprocess_time += t2 - t1
                        
                        outputs = outputs.cpu().detach().numpy()
                        transformed_batches = transformed_batches.cpu().detach().numpy()
                        lengths_list = []
                        for t in range(8):
                            lengths = evaluator.eval_batch(transformed_batches[t:t+1], outputs,  REST_model.actor.degree)
                            lengths_list.append(lengths)
                        aligned_generated_routing_twl += min(lengths_list)

                print(f'batch: {batch_id}/{num_batch}, ' + \
                        f'generated wl: {aligned_generated_routing_twl}, best wl: {aligned_real_routing_twl}, ' + \
                        f'prop: {aligned_generated_routing_twl/aligned_real_routing_twl}, ' + \
                        f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                        f'total time: {generation_time + postprocess_time}', end='\r')
                batch_id += 1
        print()
        print(f'batch: {batch_id}/{num_batch}, ' + \
                f'generated wl: {aligned_generated_routing_twl}, best wl: {aligned_real_routing_twl}, ' + \
                f'prop: {aligned_generated_routing_twl/aligned_real_routing_twl}, ' + \
                f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                f'total time: {generation_time + postprocess_time}')
            
    elif args.case == 'ISPD':
        data_config = OmegaConf.load(f"configs/test-{args.case}.yaml")
        data_module = DInterface(data_config.data, batch_size=batch_size)
        aligned_real_routing_twl = 0.
        aligned_generated_routing_twl = 0.
        generation_time = 0.
        postprocess_time = 0.
        batch_id = 0
        num_batch = len(data_module.all_dataloader())
        print('num_batch: ', num_batch)
        
        with torch.no_grad():
            for batch in data_module.all_dataloader():
                c = batch['condition'].to(device)
                x = batch['route'].to(device)
                
                real_pin = c[:, :1, :, :].detach()
                real_routing_length = torch.sum(x > 0.6)
                
                aligned_real_routing_twl += (real_routing_length.sum().item()-batch_size)
                
                t0 = time()
                if args.model == 'DDPM':
                    # model.model.ddim_timesteps = 75
                    x_gen = model.model.ddim_sample(c, guide_w=1.5)
                elif args.model == 'GAN':
                    x_gen = model(c)
                elif args.model == 'VAE':
                    z = torch.randn(c.shape[0], config.model.vae_config.params.latent_size).to(device)
                    x_gen = model.vae.inference(z, c)
                t1 = time()
    
                generation_time += t1 - t0
                x_mask = torch.sign((x_gen.sum(2).unsqueeze(2) + x_gen.sum(3).unsqueeze(-1))).to(device)
                x_mask = (x_mask + 1)/2. # to [0,1]
                x_gen = x_gen * x_mask
                
                # also add rectangular mask
                for batch_j in range(real_pin.shape[0]):
                    pin_index = torch.where(real_pin[batch_j] == 1.)
                    mask_route = torch.zeros_like(real_pin[batch_j])
                    mask_route[0, pin_index[1].min() : pin_index[1].max() + 1, 
                           pin_index[2].min() : pin_index[2].max() + 1] =  1.
                    x_gen[batch_j] = x_gen[batch_j] * mask_route
                
                x_vertex = torch.max(x_gen, real_pin)
                
                t1 = time()
                ## Geosteiner
                if data_config.data.rsmt == 'GEO':
                    for test_case in x_vertex[:, 0, :, :]:
                        tc = test_case.detach().cpu().numpy()
                        if args.model in ['DDPM', 'VAE']:
                            tc = np.array(np.where(tc >= 5/6)).T
                        elif args.model in ['GAN']:
                            tc = np.array(np.where(tc >= 1/2)).T
                        gst_length, sp_list, edge_list = evaluator.gst_rsmt(tc)
                        aligned_generated_routing_twl += gst_length
                    t2 = time()
                    postprocess_time += t2 - t1
                print(f'batch: {batch_id}/{num_batch}, ' + \
                        f'generated wl: {aligned_generated_routing_twl}, best wl: {aligned_real_routing_twl}, ' + \
                        f'prop: {aligned_generated_routing_twl/aligned_real_routing_twl}, ' + \
                        f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                        f'total time: {generation_time + postprocess_time}', end='\r')
                batch_id += 1
        print()
        print(f'batch: {batch_id}/{num_batch}, ' + \
                f'generated wl: {aligned_generated_routing_twl}, best wl: {aligned_real_routing_twl}, ' + \
                f'prop: {aligned_generated_routing_twl/aligned_real_routing_twl}, ' + \
                f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                f'total time: {generation_time + postprocess_time}')
        
if __name__ == '__main__':
    args = get_args()
    inference(args)