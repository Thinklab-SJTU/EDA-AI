import os
import sys
import numpy as np
from time import time
import torch
import torch.nn as nn
import torch.optim.lr_scheduler as lrs
import pytorch_lightning as pl
from torchvision.utils import save_image, make_grid
import matplotlib.pyplot as plt
from functools import reduce

from REST_tool.REST_utils import Evaluator ### REST
from utils import instantiate_from_config
from utils.obstacle_utils import *
from utils.rsmt_utils import *
from utils.plot_utils import *

class ACInterface(pl.LightningModule):
    def __init__(self, 
                 model_config, 
                 model_name, 
                 img_dir, 
                 batch_size, 
                 min_degree, 
                 max_degree, 
                 min_obstacle, 
                 max_obstacle, 
                 obstacle_weight, 
                 transformation, 
                 seed):
        super().__init__()
        self.save_hyperparameters()
        self.seed = str(seed)
        self.values = dict() # log_dict
        
        self.actor = instantiate_from_config(model_config.actor_config)
        self.critic = instantiate_from_config(model_config.critic_config)
        
        self.batch_size = batch_size
        self.model_name = model_name
        self.lr = model_config.lr
        self.min_degree = min_degree
        self.max_degree = max_degree
        self.min_obstacle = min_obstacle
        self.max_obstacle = max_obstacle
        self.obstacle_weight = obstacle_weight
        self.img_dir = img_dir
        self.evaluator = Evaluator(path=os.getcwd() + '/REST_tool/algorithms/libeval.so') # We borrow the evaluator from REST
        self.mse_loss = torch.nn.MSELoss()
        self.transformation = transformation

    def forward(self, input_batch, input_obstacle_points_batch, trans=False):
        outputs, log_probs = self.actor(input_batch, input_obstacle_points_batch, self.max_obstacle, trans)
        return outputs, log_probs

    def training_step(self, batch, batch_idx):
        self.actor.train()
        self.critic.train()
        input_batch = batch[0].to(self.device)
        input_obstacles_batch = batch[1].to(self.device)
        input_obstacle_points_batch = batch[2].to(self.device)
        outputs, log_probs = self(input_batch, input_obstacle_points_batch)
        predictions_oa_length, predictions_ol = self.critic(input_batch, input_obstacle_points_batch, self.max_obstacle)
        lengths_oa, overlaps, _, _ = self.evaluator.eval_obstacle_batch(input_batch.cpu().detach().numpy(), 
                input_obstacles_batch.cpu().detach().numpy(), outputs.cpu().detach().numpy(), self.max_degree)
        length_oa_tensor = torch.tensor(lengths_oa, dtype=torch.float, device=self.device)
        overlap_tensor = torch.tensor(overlaps, dtype=torch.float, device=self.device)

        with torch.no_grad():
            disadvantage_oa = length_oa_tensor - predictions_oa_length + (overlap_tensor) * self.obstacle_weight
        actor_oa_loss = torch.mean(disadvantage_oa * log_probs)
        critic_oa_length_loss = self.mse_loss(predictions_oa_length, length_oa_tensor)
        if self.max_obstacle > 0:
            critic_ob_loss = self.mse_loss(predictions_ol, overlap_tensor)
            loss = actor_oa_loss + (critic_oa_length_loss + critic_ob_loss)
        else:
            loss = actor_oa_loss + critic_oa_length_loss
        self.log_util(loss, 'train_loss')

        return loss
    
    def validation_step(self, batch, batch_idx):
        self.actor.eval()
        self.critic.eval()
        input_batch = batch[0].to(self.device)
        input_obstacles_batch = batch[1].to(self.device)
        input_obstacle_points_batch = batch[2].to(self.device)
        with torch.no_grad():
            outputs, log_probs = self(input_batch, input_obstacle_points_batch, True)
            predictions_length, predictions_ol = self.critic(input_batch, input_obstacle_points_batch, self.max_obstacle)
        lengths, overlaps, _, _  = self.evaluator.eval_obstacle_batch(input_batch.cpu().detach().numpy(), 
                input_obstacles_batch.cpu().detach().numpy(), outputs.cpu().detach().numpy(), self.max_degree)
        length_oa_tensor = torch.tensor(lengths, dtype=torch.float, device=self.device)
        overlap_tensor = torch.tensor(overlaps, dtype=torch.float, device=self.device)

        with torch.no_grad():
            disadvantage_oa = length_oa_tensor - predictions_length + (overlap_tensor) * self.obstacle_weight
        actor_oa_loss = torch.mean(disadvantage_oa * log_probs)
        critic_oa_length_loss = self.mse_loss(predictions_length, length_oa_tensor)
        if self.max_obstacle > 0:
            critic_ob_loss = self.mse_loss(predictions_ol, overlap_tensor)
            loss = actor_oa_loss + (critic_oa_length_loss + critic_ob_loss)
        else:
            loss = actor_oa_loss + critic_oa_length_loss
        self.log_util(loss, 'val_loss')
        self.log_util(lengths.mean(), 'val_lengths')
        self.log_util(overlaps.mean(), 'val_overlaps')
        self.log_util(lengths.mean() + overlaps.mean() * self.obstacle_weight, 'val_metric')
        
        # draw results
        if batch_idx == 0:
            num_samples = 10
            all_lengths = []
            all_overlaps = []
            all_new_inputs = []
            all_new_outputs = []
            all_new_obstacles = []
            all_t = []
            best_lengths = [1e9 for i in range(num_samples)]
            best_overlaps = [1e9 for i in range(num_samples)]
            best_t = [0 for i in range(num_samples)]
            best_new_inputs = [[] for i in range(num_samples)]
            best_new_outputs = [[] for i in range(num_samples)]
            best_new_obstacles = [[] for i in range(num_samples)]
            
            for t in range(self.transformation):
                transformed_batch = transform_inputs(input_batch, t)
                transformed_obstacle_batch = transform_obstacles(input_obstacles_batch, t)
                transformed_obstacle_points_batch = obstacle_to_points(transformed_obstacle_batch)
                with torch.no_grad():
                    outputs, _ = self.actor(transformed_batch, transformed_obstacle_points_batch, self.max_obstacle, True)
                outputs = outputs.cpu().detach().numpy()
                lengths, overlaps, new_inputs, new_outputs = self.evaluator.eval_obstacle_batch(transformed_batch.detach().cpu().numpy(), transformed_obstacle_batch.detach().cpu().numpy(), outputs, self.max_degree)
                for i in range(num_samples):
                    if lengths[i] + overlaps[i] * 10e6 < best_lengths[i] + best_overlaps[i] * 10e6:
                        best_lengths[i] = lengths[i]
                        best_overlaps[i] = overlaps[i]
                        best_new_inputs[i] = new_inputs[i]
                        best_new_outputs[i] = new_outputs[i]
                        best_new_obstacles[i] = transformed_obstacle_batch[i].detach().cpu().numpy()
                        best_t[i] = t
                    
            all_lengths.append(best_lengths)
            all_overlaps.append(best_overlaps)
            all_new_outputs.append(best_new_outputs)
            all_new_inputs.append(best_new_inputs)
            all_new_obstacles.append(best_new_obstacles)
            all_t.append(best_t)

            all_lengths = np.concatenate(all_lengths, 0)
            all_overlaps = np.concatenate(all_overlaps, 0)
            all_new_outputs = reduce(lambda x,y: x+y, all_new_outputs) 
            all_new_inputs = reduce(lambda x,y: x+y, all_new_inputs) 
            all_new_obstacles = np.concatenate(all_new_obstacles, 0)
            all_t = np.concatenate(all_t, 0)

            self.draw_results(all_lengths, all_overlaps, all_new_outputs, all_new_inputs, all_new_obstacles)

        return lengths, overlaps

    def test_step(self, batch, batch_idx):
        self.actor.eval()
        self.critic.eval()
        input_batch = batch[0].to(self.device)
        input_obstacles_batch = batch[1].to(self.device)

        # visualize results
        if batch_idx == 0:
            num_samples = 50
            print(f"Visualize {num_samples} samples.")
            all_lengths = []
            all_overlaps = []
            all_new_inputs = []
            all_new_outputs = []
            all_new_obstacles = []
            all_t = []
            best_lengths = np.ones(num_samples, dtype=float) * 1e9
            best_overlaps = np.ones(num_samples, dtype=float) * 1e9
            best_t = np.zeros(num_samples, dtype=float)
            best_new_inputs = [[]] * num_samples
            best_new_outputs = [[]] * num_samples
            best_new_obstacles = [[]] * num_samples
            test_case = input_batch[:num_samples]
            for t in range(self.transformation):
                # if t != 2:
                #     continue
                transformed_batch = transform_inputs(test_case, t)
                transformed_obstacle_batch = transform_obstacles(input_obstacles_batch[:num_samples], t)
                transformed_obstacle_points_batch = obstacle_to_all_points(transformed_obstacle_batch)
                with torch.no_grad():
                    outputs, _ = self.actor(transformed_batch, transformed_obstacle_points_batch, self.max_obstacle, True)
                outputs = outputs.cpu().detach().numpy()
                # print("transformed_batch.detach().cpu().numpy()[0]: ", transformed_batch.detach().cpu().numpy()[0])
                # print("transformed_obstacle_batch.detach().cpu().numpy()[0]: ", transformed_obstacle_batch.detach().cpu().numpy()[0])
                # print("outputs: ", outputs[0])
                # sys.exit()
                lengths, overlaps, new_inputs, new_outputs = self.evaluator.eval_obstacle_batch(transformed_batch.detach().cpu().numpy(), transformed_obstacle_batch.detach().cpu().numpy(), outputs, self.max_degree)
                for i in range(num_samples):
                    if lengths[i] + overlaps[i] * 1e6 < best_lengths[i] + best_overlaps[i] * 1e6:
                        best_lengths[i] = lengths[i]
                        best_overlaps[i] = overlaps[i]
                        best_new_inputs[i] = new_inputs[i]
                        best_new_outputs[i] = new_outputs[i]
                        best_new_obstacles[i] = transformed_obstacle_batch[i].detach().cpu().numpy()
                        best_t[i] = t
            # print("best_t[0]: ", best_t[0])
            # sys.exit()
            all_lengths.append(best_lengths)
            all_overlaps.append(best_overlaps)
            all_new_outputs.append(best_new_outputs)
            all_new_inputs.append(best_new_inputs)
            all_new_obstacles.append(best_new_obstacles)
            all_t.append(best_t)

            all_lengths = np.concatenate(all_lengths, 0)
            all_overlaps = np.concatenate(all_overlaps, 0)
            all_new_outputs = reduce(lambda x,y: x+y, all_new_outputs) 
            all_new_inputs = reduce(lambda x,y: x+y, all_new_inputs) 
            all_new_obstacles = np.concatenate(all_new_obstacles, 0)
            all_t = np.concatenate(all_t, 0)

            self.draw_results(all_lengths, all_overlaps, all_new_outputs, all_new_inputs, all_new_obstacles, test=True)
            print(f"Finished.")
            # sys.exit()

        batch_size = input_batch.shape[0]
        all_lengths = []
        all_overlaps = []
        best_lengths = np.ones(batch_size, dtype=float) * 1e9
        best_overlaps = np.ones(batch_size, dtype=float) * 1e9
        best_t = np.zeros(batch_size, dtype=int)
        for t in range(self.transformation):
            transformed_batch = transform_inputs(input_batch, t)
            transformed_obstacle_batch = transform_obstacles(input_obstacles_batch, t)
            transformed_obstacle_points_batch = obstacle_to_all_points(transformed_obstacle_batch)
            start_time = time()
            with torch.no_grad():
                outputs, _ = self.actor(transformed_batch, transformed_obstacle_points_batch, self.max_obstacle, True)
            oarest_time = time() - start_time
            self.total_oarest_time += oarest_time
            outputs = outputs.cpu().detach().numpy()
            if self.max_obstacle == 0 and self.min_degree == self.max_degree:
                lengths = self.evaluator.eval_batch(transformed_batch.detach().cpu().numpy(), 
                                                    outputs, self.max_degree)
                overlaps = lengths.copy() * 0.0
            else:
                lengths, overlaps, new_inputs, new_outputs = self.evaluator.eval_obstacle_batch(transformed_batch.detach().cpu().numpy(), 
                                                                                                transformed_obstacle_batch.detach().cpu().numpy(), outputs, self.max_degree)
            for i in range(batch_size):
                if lengths[i] + overlaps[i] * 1e6 < best_lengths[i] + best_overlaps[i] * 1e6:
                    best_lengths[i] = lengths[i]
                    best_overlaps[i] = overlaps[i]
                    best_t[i] = t
        all_lengths.append(best_lengths)
        all_overlaps.append(best_overlaps)
        all_lengths = np.concatenate(all_lengths, 0)
        all_overlaps = np.concatenate(all_overlaps, 0)
        self.oarest_lengths_accum.append(all_lengths.copy())
        self.oarest_overlaps_accum.append(all_overlaps.copy())

        # GeoSt
        start_time = time()
        gst_lengths = []
        gst_overlaps = []
        inputs = input_batch.detach().cpu().numpy()
        obstacles = input_obstacles_batch.detach().cpu().numpy()
        for i in range(len(inputs)):
            test_case = inputs[i]
            test_obstacle = obstacles[i]
            test_case = test_case[test_case[:, 0] != -1]
            # gst_length, _, _ = self.evaluator.gst_rsmt(test_case)
            gst_length, gst_overlap, _, _ = self.evaluator.gst_obstacle_rsmt(test_case, test_obstacle)
            gst_lengths.append(gst_length)
            gst_overlaps.append(gst_overlap)
        gst_lengths = np.array(gst_lengths)
        gst_overlaps = np.array(gst_overlaps)
        rsmt_time = time() - start_time
        self.total_rsmt_time += rsmt_time
        self.gst_lengths_accum.append(gst_lengths)
        self.gst_overlaps_accum.append(gst_overlaps)

    def log_util(self, loss, name='loss'):
        self.values[name] = loss
        self.log_dict(self.values, logger=True, prog_bar=True, on_step=True, on_epoch=True, 
                      batch_size=self.batch_size)
        
    def on_train_epoch_end(self):
        torch.nn.utils.clip_grad_norm_(self.actor.parameters(), 1.)
        torch.nn.utils.clip_grad_norm_(self.critic.parameters(), 1.)

    def on_test_epoch_start(self):
        self.total_oarest_time = 0
        self.total_rsmt_time = 0
        self.oarest_lengths_accum = []
        self.oarest_overlaps_accum = []
        self.gst_lengths_accum = [] 
        self.gst_overlaps_accum = []
        self.total_samples = 0

    def on_test_epoch_end(self):
        self.log("test_oarest_time", self.total_oarest_time, prog_bar=True)
        all_oarest_lengths = np.concatenate(self.oarest_lengths_accum, axis=0)
        all_oarest_overlaps = np.concatenate(self.oarest_overlaps_accum, axis=0)
        oarest_mean = all_oarest_lengths.mean()
        overlap_mean = all_oarest_overlaps.mean()
        self.log("test_oarest_length", oarest_mean)
        self.log("test_oarest_overlap", overlap_mean)
        if len(self.gst_lengths_accum) != 0:
            self.log("test_rsmt_time", self.total_rsmt_time, prog_bar=True)
            all_gst_lengths = np.concatenate(self.gst_lengths_accum, axis=0)
            all_gst_overlaps = np.concatenate(self.gst_overlaps_accum, axis=0)
            gst_mean = all_gst_lengths.mean()
            gst_overlap_mean = all_gst_overlaps.mean()
            self.log("test_geost_length", gst_mean)
            self.log("test_geost_overlap", gst_overlap_mean)
            self.log("test_length_error", ((all_oarest_lengths-all_gst_lengths)/all_gst_lengths).mean())
        self.log("num_samples", len(all_oarest_lengths))

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(list(self.actor.parameters()) + list(self.critic.parameters()), lr=self.lr, eps=1e-5)
        
        return [optimizer], []

    def draw_results(self, all_lengths, all_overlaps, all_new_outputs, all_new_inputs, all_new_obstacles, test=False):
        if test:
            dir_name = f"{self.img_dir}test_{self.model_name}_degree-[{self.min_degree},{self.max_degree}]_obstacle-[{self.min_obstacle},{self.max_obstacle}]_seed-{self.seed}/"
        else:
            dir_name = f"{self.img_dir}{self.model_name}_degree-[{self.min_degree},{self.max_degree}]_obstacle-[{self.min_obstacle},{self.max_obstacle}]_seed-{self.seed}/"
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        plt.rcParams['font.family'] = 'DejaVu Serif'
        for i in range(len(all_lengths)):
            fig = plt.figure(figsize=(10, 4.6))
            ax1 = plt.subplot(1, 2, 1)
            # Optimal RSMT
            # print(all_new_inputs[i], all_new_outputs[i], all_new_obstacles[i])
            gst_length, sps, edges = self.evaluator.gst_rsmt(all_new_inputs[i][:self.max_degree])
            ax1 = plot_obstacle_gst_rsmt(all_new_inputs[i][:self.max_degree], sps, edges, all_new_obstacles[i], ax1)
            if len(sps) > 0:
                overlap = count_overlap(np.concatenate([all_new_inputs[i][:self.max_degree], sps], axis=0), all_new_obstacles[i], reduce(lambda x,y: x+y, edges))
            else:
                overlap = count_overlap(all_new_inputs[i][:self.max_degree], all_new_obstacles[i], reduce(lambda x,y: x+y, edges))
            plt.annotate('Optimal Length: ' + str(round(gst_length, 3)) + ', Overlap: ' + str(overlap), (-0.04, -0.04))

            ax2 = plt.subplot(1, 2, 2)
            # OAREST solution
            ax2 = plot_obstacle_rest(all_new_inputs[i], all_new_outputs[i], all_new_obstacles[i], ax2, self.max_degree)
            plt.annotate('OAREST Length: ' + str(round(all_lengths[i], 3)) + ', Overlap: ' + str(int(all_overlaps[i])), (-0.04, -0.04))

            fig.savefig(dir_name + \
                        f"epoch-{self.current_epoch}_sample-{i}_gst-{round(gst_length, 3)}_oarest-{round(all_lengths[i], 3)}.pdf", bbox_inches='tight')