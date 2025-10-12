import sys
import numpy as np
from time import time
import torch
import torch.nn as nn
import torch.nn.functional as F
import pytorch_lightning as pl

from torch.distributions.categorical import Categorical
from .template_nets import Embedder, Pointer, Glimpse
from .self_attn import Encoder

def check_orthogonal_path_intersects_batch(start_points, end_points, obstacles):
    """
    start_points: [batch_size, 2] start points
    end_points: [batch_size, 2] end points
    obstacles: [batch_size, num_obstacles, 4, 2] obstacle vertices
    returns: [batch_size, num_obstacles] whether overlapping
    """
    batch_size, num_obstacles = obstacles.shape[:2]
    
    # AABB
    path_min_x = torch.minimum(start_points[:, 0], end_points[:, 0]).unsqueeze(1)
    path_max_x = torch.maximum(start_points[:, 0], end_points[:, 0]).unsqueeze(1)
    path_min_y = torch.minimum(start_points[:, 1], end_points[:, 1]).unsqueeze(1)
    path_max_y = torch.maximum(start_points[:, 1], end_points[:, 1]).unsqueeze(1)
    
    obs_min_x = obstacles[..., 0].min(dim=2)[0]
    obs_max_x = obstacles[..., 0].max(dim=2)[0]
    obs_min_y = obstacles[..., 1].min(dim=2)[0]
    obs_max_y = obstacles[..., 1].max(dim=2)[0]
    
    aabb_overlap = (path_min_x < obs_max_x) & (path_max_x > obs_min_x) & \
                  (path_min_y < obs_max_y) & (path_max_y > obs_min_y)
    
    result = torch.zeros((batch_size, num_obstacles), dtype=torch.bool, device=start_points.device)
    
    check_indices = torch.nonzero(aabb_overlap, as_tuple=True)
    if len(check_indices[0]) == 0:
        return result
    
    def check_line_segment(p1, p2, obstacles_to_check):
        p1 = p1[check_indices[0]]
        p2 = p2[check_indices[0]]
        obs = obstacles[check_indices[0], check_indices[1]]
        
        seg_dir = p2 - p1
        intersects = torch.zeros(len(check_indices[0]), dtype=torch.bool, device=p1.device)
        
        for i in range(4):
            edge_start = obs[:, i]
            edge_end = obs[:, (i + 1) % 4]
            
            edge_dir = edge_end - edge_start
            
            # a×b = a.x*b.y - a.y*b.x
            cross_1 = seg_dir[:, 0] * (edge_start[:, 1] - p1[:, 1]) - seg_dir[:, 1] * (edge_start[:, 0] - p1[:, 0])
            cross_2 = seg_dir[:, 0] * (edge_end[:, 1] - p1[:, 1]) - seg_dir[:, 1] * (edge_end[:, 0] - p1[:, 0])
            cross_3 = edge_dir[:, 0] * (p1[:, 1] - edge_start[:, 1]) - edge_dir[:, 1] * (p1[:, 0] - edge_start[:, 0])
            cross_4 = edge_dir[:, 0] * (p2[:, 1] - edge_start[:, 1]) - edge_dir[:, 1] * (p2[:, 0] - edge_start[:, 0])
            
            intersects |= (cross_1 * cross_2 < 0) & (cross_3 * cross_4 < 0)
            
        return intersects
    
    mid_points = torch.stack([start_points[:, 0], end_points[:, 1]], dim=1)
    intersects_1 = check_line_segment(start_points, mid_points, obstacles)
    intersects_2 = check_line_segment(mid_points, end_points, obstacles)
    
    result[check_indices[0], check_indices[1]] = intersects_1 | intersects_2
    
    return result

def find_nearest_vertices_batch(points, obstacles, mask):
    """
    points: [batch_size, 2]
    obstacles: [batch_size, num_obstacles, 4, 2]
    returns: vertex_indices [batch_size, num_obstacles], distances [batch_size, num_obstacles]
    """
    batch_size, num_obstacles = obstacles.shape[:2]
    mask_expanded = mask.view(batch_size, num_obstacles, 4, 1).expand_as(obstacles)
    
    obstacles = torch.where(mask_expanded == 1, torch.full_like(obstacles, float('inf')), obstacles)
    points = points.unsqueeze(1).unsqueeze(1).expand(-1, num_obstacles, 4, -1)
    
    distances = torch.sum(torch.abs(obstacles - points), dim=-1)  # [batch_size, num_obstacles, 4]
    min_distances, vertex_indices = distances.min(dim=-1)  # [batch_size, num_obstacles]
    return vertex_indices, min_distances

class Actor(pl.LightningModule):
    def __init__(self, max_degree):
        super().__init__()
        self.save_hyperparameters()
        self.max_degree = max_degree

        # embedder args
        self.d_input = 2
        self.d_model = 128
        self.embedder = Embedder(max_degree, self.d_input, self.d_model)
        # self.embedder_ob_points = Embedder(max_degree, self.d_input, self.d_model)
        # self.embedder_ob = nn.Linear(3*self.d_model, self.d_model, bias = False)

        # encoder args
        self.num_stacks = 3
        self.num_heads = 16
        self.d_k = 16
        self.d_v = 16
        # feedforward layer inner
        self.d_inner = 512

        self.encoder = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        # self.encoder_ob_points = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        # self.encoder_ob = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        # self.encoder_multi = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)

        # decoder args
        self.d_unit = 256
        self.d_query = 360
        self.conv1d_r = nn.Conv1d(self.d_model, self.d_unit, 1)
        self.conv1d_x = nn.Conv1d(self.d_model, self.d_unit, 1)
        self.conv1d_y = nn.Conv1d(self.d_model, self.d_unit, 1)

        self.start_ptr = Pointer(self.d_query, self.d_unit)
        self.q_l1 = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_l2 = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_l3 = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_lx = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_ly = nn.Linear(self.d_model, self.d_query, bias=False)
        self.relu = nn.ReLU()
        self.ctx_linear = nn.Linear(self.d_query, self.d_query, bias=False)

        self.ptr1 = Pointer(self.d_query, self.d_unit)
        self.ptr2 = Pointer(self.d_query, self.d_unit)

        # self.to(device)
        self.train()

    def forward(self, inputs, input_obstacle_points_batch=None, max_obstacle=0, deterministic=False):
        if max_obstacle == 0:
            return self.forward_without_obstacle(inputs, deterministic)
        else:
            return self.forward_with_obstacle(inputs, input_obstacle_points_batch, deterministic)
    
    def forward_without_obstacle(self, inputs, deterministic=False):
        inputs_tensor = torch.as_tensor(inputs).to(dtype=torch.float, device=self.device)
        # create mask for inputs_tensor
        with torch.no_grad():
            inputs_mask = (inputs_tensor[:, :, 0] == -1).float() # only check the first axis
            inputs_mask_expand = inputs_mask.unsqueeze(-1)
        batch_size = inputs_tensor.size(0)
        degree_range = torch.arange(self.max_degree, device=self.device)
        torch.cuda.synchronize()
        embedings = self.embedder(inputs_tensor).masked_fill(inputs_mask_expand == 1, 0) # mask the -1 values
        encodings = self.encoder(embedings, None)
        
        encodings_range = torch.arange(encodings.size(0))
        enc_r = self.conv1d_r(encodings.transpose(1, 2)).transpose(1, 2)
        enc_x = self.conv1d_x(encodings.transpose(1, 2)).transpose(1, 2)
        enc_y = self.conv1d_y(encodings.transpose(1, 2)).transpose(1, 2)
        enc_xy = torch.cat([enc_x, enc_y], dim=1)
        # create indicies
        batch_indices = torch.arange(batch_size, device=self.device)

        visited = torch.zeros([batch_size, self.max_degree], dtype=torch.bool, device=self.device)

        # find the positions of the first -1
        first_neg_mask = (inputs_mask.cumsum(dim=1) == 1).float()
        # first_neg_positions = first_neg_mask.argmax(dim=1)
        row_has_neg = (first_neg_mask.sum(dim=1) > 0)
        first_neg_positions = first_neg_mask.argmax(dim=1)
        first_neg_positions = first_neg_positions.masked_fill(~row_has_neg, self.max_degree)
        
        valid_inputs = first_neg_positions < inputs_tensor.size(1)  # find valid positions
        if valid_inputs.any():
            # create range mask
            positions = first_neg_positions[valid_inputs]
            row_idx = batch_indices[valid_inputs]
            # fill visited
            visited[row_idx] |= degree_range >= positions.unsqueeze(1)
        # create mask for obstacles
        init_mask = visited.detach().clone().to(self.device)

        indexes, log_probs = [], []

        start_logits = self.start_ptr(enc_r,
            torch.zeros([batch_size, self.d_query], dtype=torch.float).to(self.device), visited)
        distr = Categorical(logits=start_logits)
        if deterministic:
            start_idx = start_logits.max(dim=-1)[1]
        else:
            start_idx = distr.sample()
        visited_indices = start_idx.unsqueeze(-1)
        visited.scatter_(1, visited_indices, True)
        log_probs.append(distr.log_prob(start_idx))

        q1 = q2 = qx = qy = encodings[torch.arange(batch_size), start_idx]
        context = torch.zeros([batch_size, self.d_query], device=self.device)

        for step in range(self.max_degree - 1):
            residual = self.q_l1(q1).add_(self.q_l2(q2)).add_(self.q_lx(qx)).add_(self.q_ly(qy))
            context = torch.max(context, self.ctx_linear(self.relu(residual)))
            first_q = residual.add_(context)
            first_query = self.relu(first_q)
            # first_idx
            logits = self.ptr1(enc_r, first_query, visited | init_mask)

            b = logits.size(0)
            all_neg_inf = (logits == float('-inf')).all(dim=1).to(self.device)
            first_idx = torch.full((b,), 0, dtype=torch.long).to(self.device)

            log_prob = torch.zeros(b).to(self.device)

            valid_mask = ~all_neg_inf
            if valid_mask.any():
                valid_logits = logits[valid_mask]
                distr = Categorical(logits=valid_logits)
                if deterministic:
                    samples = valid_logits.max(dim=-1)[1]
                else:
                    samples = distr.sample()  
                first_idx[valid_mask] = samples
                log_prob[valid_mask] = distr.log_prob(samples)

            log_probs.append(log_prob)

            # second_idx
            q3 = encodings[encodings_range, first_idx]
            second_query = self.relu(first_q.add_(self.q_l3(q3)))
            unvisited = (~visited) | init_mask
            unvisited = unvisited.repeat(1, 2)
            logits = self.ptr2(enc_xy, second_query, unvisited)
            b = logits.size(0)
            all_neg_inf = (logits == float('-inf')).all(dim=1).to(self.device)
            idxs = torch.full((b,), 1, dtype=torch.long).to(self.device)

            log_prob = torch.zeros(b).to(self.device)

            valid_mask = ~all_neg_inf
            if valid_mask.any():
                valid_logits = logits[valid_mask]
                distr = Categorical(logits=valid_logits)
                if deterministic:
                    samples = valid_logits.max(dim=-1)[1]
                else:
                    samples = distr.sample()
                idxs[valid_mask] = samples
                log_prob[valid_mask] = distr.log_prob(samples)
            log_probs.append(log_prob)
            second_idx = idxs % self.max_degree
            
            sec_dir = torch.div(idxs, self.max_degree, rounding_mode='floor')
            fir_dir = 1 - sec_dir
            with torch.no_grad():
                x_idx = first_idx * sec_dir + second_idx * fir_dir
                y_idx = first_idx * fir_dir + second_idx * sec_dir
            indexes.append(x_idx)
            indexes.append(y_idx)

            # update visited
            visited.scatter_(1, first_idx.unsqueeze(-1), True)

            # update query
            q1 = q3
            q2 = encodings[encodings_range, second_idx]
            qx = encodings[encodings_range, x_idx]
            qy = encodings[encodings_range, y_idx]
        indexes = torch.stack(indexes, -1)
        log_probs = sum(log_probs)
        return indexes, log_probs

    def forward_with_obstacle(self, inputs, input_obstacle_points_batch, deterministic=False):
        inputs_tensor = torch.as_tensor(inputs).to(dtype=torch.float, device=self.device)
        input_obstacle_points_tensor = torch.as_tensor(input_obstacle_points_batch).to(dtype=torch.float, device=self.device)
        inputs_with_obstacle_tensor = torch.cat([inputs_tensor, input_obstacle_points_tensor], dim=1)

        with torch.no_grad():
            inputs_mask = (inputs_tensor[:, :, 0] == -1) # only check the first axis
            obstacle_mask = (input_obstacle_points_tensor[:, :, 0] == -1)
            invalid_mask = torch.cat([inputs_mask, obstacle_mask], dim=1).float()
            invalid_mask_expand = invalid_mask.unsqueeze(-1)
        batch_size = inputs_tensor.size(0)
        self.total_length = self.max_degree + input_obstacle_points_tensor.shape[1]
        pos_range = torch.arange(self.total_length, device=inputs_tensor.device).unsqueeze(0)
        degree_range = torch.arange(self.max_degree, device=self.device)
        torch.cuda.synchronize()
        embedings = self.embedder(inputs_with_obstacle_tensor).masked_fill(invalid_mask_expand == 1, 0) # mask the -1 values
        encodings = self.encoder(embedings, None)
        encodings_range = torch.arange(encodings.size(0))
        enc_r = self.conv1d_r(encodings.transpose(1, 2)).transpose(1, 2)
        enc_x = self.conv1d_x(encodings.transpose(1, 2)).transpose(1, 2)
        enc_y = self.conv1d_y(encodings.transpose(1, 2)).transpose(1, 2)
        enc_xy = torch.cat([enc_x, enc_y], dim=1)

        # create indicies
        batch_indices = torch.arange(batch_size, device=self.device)

        self.visited = torch.zeros([batch_size, self.total_length], dtype=torch.bool, device=self.device)
        start_visited = torch.zeros([batch_size, self.total_length], dtype=torch.bool, device=self.device)
        self.init_mask = torch.zeros([batch_size, self.total_length], dtype=torch.bool, device=self.device)
        start_visited[:, self.max_degree:] = 1
        self.unactivated_obstacle = torch.zeros([batch_size, self.total_length], dtype=torch.bool, device=self.device)
        self.unactivated_obstacle[:, self.max_degree:] = 1
        self.init_unactivated_obstacle = self.unactivated_obstacle.detach().clone()

        first_neg_mask = (inputs_mask.float().cumsum(dim=1) == 1).float()
        row_has_neg = (first_neg_mask.sum(dim=1) > 0)
        first_neg_positions = first_neg_mask.argmax(dim=1)
        first_neg_positions = first_neg_positions.masked_fill(~row_has_neg, self.max_degree)
        valid_inputs = first_neg_positions < inputs_tensor.size(1)  # find valid positions
        if valid_inputs.any():
            # create range mask
            positions = first_neg_positions[valid_inputs]
            row_idx = batch_indices[valid_inputs]
            # fill start_visited
            start_visited[row_idx] |= pos_range >= positions.unsqueeze(1)
            # fill visited (until max_degree)
            self.init_mask[row_idx, :self.max_degree] |= degree_range >= positions.unsqueeze(1)
        
        # mask all obstacles at first
        first_neg_mask_obstacle = (obstacle_mask.float().cumsum(dim=1) == 1).float()
        row_has_neg = (first_neg_mask_obstacle.sum(dim=1) > 0)
        first_neg_positions_obstacle = first_neg_mask_obstacle.argmax(dim=1)
        first_neg_positions_obstacle = first_neg_positions_obstacle.masked_fill(~row_has_neg, self.total_length)
        valid_obstacles = first_neg_positions_obstacle < input_obstacle_points_tensor.size(1)
        if valid_obstacles.any():
            positions = first_neg_positions_obstacle[valid_obstacles]
            row_idx = batch_indices[valid_obstacles]
            obstacle_start = self.max_degree + positions
            self.init_mask[row_idx] |= pos_range >= obstacle_start.unsqueeze(1)

        indexes, log_probs = [], []
        start_logits = self.start_ptr(enc_r,
            torch.zeros([batch_size, self.d_query], dtype=torch.float).to(self.device), start_visited)
        distr = Categorical(logits=start_logits)
        if deterministic:
            start_idx = start_logits.max(dim=-1)[1]
        else:
            start_idx = distr.sample()
        self.visited.scatter_(1, start_idx.unsqueeze(-1), True)
        log_probs.append(distr.log_prob(start_idx))

        q1 = q2 = qx = qy = encodings[torch.arange(batch_size), start_idx]
        context = torch.zeros([batch_size, self.d_query], device=self.device)
        for step in range(self.total_length - 1):
            self.step = step
            residual = self.q_l1(q1).add_(self.q_l2(q2)).add_(self.q_lx(qx)).add_(self.q_ly(qy))
            context = torch.max(context, self.ctx_linear(self.relu(residual)))
            first_q = residual.add_(context)
            first_query = self.relu(first_q)

            # first_idx: select unvisited points
            # example: 
            # visited: [0, 0, 1, 0, | 0, 0, 0, 0]
            # self.init_mask: [0, 0, 0, 1, | 0, 0, 0, 0]
            # self.unactivated_obstacle: [0, 0, 0, 0, | 1, 1, 0, 1]
            # comprehensive_mask: [0, 0, 1, 1, | 1, 1, 0, 1]
            comprehensive_mask = self.visited | self.init_mask | self.unactivated_obstacle
            logits = self.ptr1(enc_r, first_query, comprehensive_mask) # value 1 means invalid
            b = logits.size(0)
            all_neg_inf = (logits == float('-inf')).all(dim=1).to(self.device)
            first_idx = torch.full((b,), self.max_degree, dtype=torch.long).to(self.device)

            log_prob = torch.zeros(b).to(self.device)
            
            valid_mask = ~all_neg_inf
            if valid_mask.any():
                valid_logits = logits[valid_mask]
                distr = Categorical(logits=valid_logits)
                if deterministic:
                    samples = valid_logits.max(dim=-1)[1]
                else:
                    samples = distr.sample()  
                first_idx[valid_mask] = samples
                log_prob[valid_mask] = distr.log_prob(samples)
            else:
                # no valid positions
                break

            log_probs.append(log_prob)
            # second_idx
            q3 = encodings[encodings_range, first_idx]
            second_query = self.relu(first_q.add_(self.q_l3(q3)))
            
            unvisited = (~self.visited)
            unvisited = unvisited.repeat(1, 2)
            logits = self.ptr2(enc_xy, second_query, unvisited)
            b = logits.size(0)
            all_neg_inf = (logits == float('-inf')).all(dim=1).to(self.device)
            idxs = torch.full((b,), self.max_degree, dtype=torch.long).to(self.device)

            log_prob = torch.zeros(b).to(self.device)

            valid_mask = ~all_neg_inf
            if valid_mask.any():
                valid_logits = logits[valid_mask]
                distr = Categorical(logits=valid_logits)
                if deterministic:
                    samples = valid_logits.max(dim=-1)[1]
                else:
                    samples = distr.sample()
                idxs[valid_mask] = samples
                log_prob[valid_mask] = distr.log_prob(samples)
            log_probs.append(log_prob)

            second_idx = idxs % self.total_length
            # update new idx if overlaps
            new_first_idx, new_second_idx, new_idxs, x_idx, y_idx, new_q3 = self.get_adjusted_path_batch(
                first_idx, second_idx, idxs, 
                inputs_with_obstacle_tensor, input_obstacle_points_tensor,
                encodings, enc_r, enc_xy, first_q, deterministic
            )
            sec_dir = torch.div(new_idxs, self.total_length, rounding_mode='floor')
            fir_dir = 1 - sec_dir
            with torch.no_grad():
                x_idx = new_first_idx * sec_dir + new_second_idx * fir_dir
                y_idx = new_first_idx * fir_dir + new_second_idx * sec_dir

            indexes.append(x_idx)
            indexes.append(y_idx)
            # update visited
            self.visited.scatter_(1, new_first_idx.unsqueeze(-1), True)

            # update query
            q1 = new_q3
            q2 = encodings[encodings_range, new_second_idx]
            qx = encodings[encodings_range, x_idx]
            qy = encodings[encodings_range, y_idx]
        indexes = torch.stack(indexes, -1)
        log_probs = sum(log_probs)
        return indexes, log_probs

    def get_adjusted_path_batch(self, first_idx, second_idx, idxs, inputs_with_obstacle_tensor, input_obstacle_points_tensor, encodings, enc_r, enc_xy, first_q, deterministic):
        """
        check and adjust path in batch
        """
        t0 = time()
        batch_size = inputs_with_obstacle_tensor.size(0)
        device = inputs_with_obstacle_tensor.device
        total_length = self.max_degree + input_obstacle_points_tensor.shape[1]
        new_idxs = idxs.clone()
        
        batch_range = torch.arange(batch_size, device=device)
        encodings_range = torch.arange(encodings.size(0))
        sec_dir = torch.div(idxs, self.total_length, rounding_mode='floor')
        fir_dir = 1 - sec_dir

        with torch.no_grad():
            x_idx = first_idx * sec_dir + second_idx * fir_dir
            y_idx = first_idx * fir_dir + second_idx * sec_dir

        first_points = inputs_with_obstacle_tensor[batch_range, x_idx]
        second_points = inputs_with_obstacle_tensor[batch_range, y_idx]
        
        num_obstacles = input_obstacle_points_tensor.size(1) // 4
        obstacles = input_obstacle_points_tensor.view(batch_size, num_obstacles, 4, 2)

        # initialization
        new_first_idx = first_idx.clone()
        new_q3 = encodings[encodings_range, first_idx]
        max_iterations = 3
        for _ in range(max_iterations):
            # check whether paths overlap with obstacles
            intersections = check_orthogonal_path_intersects_batch(
                first_points,
                second_points,
                obstacles
            )  # [batch_size, num_obstacles]
            has_intersection = intersections.any(dim=1)
            if not has_intersection.any():
                break
            
            self.unactivated_obstacle[:, self.max_degree:] = ~intersections.unsqueeze(-1).repeat(1,1,4).view(-1, intersections.shape[1]*4)

            first_query = self.relu(first_q[has_intersection])

            comprehensive_mask = self.visited[has_intersection] | self.init_mask[has_intersection] | self.unactivated_obstacle[has_intersection]

            logits = self.ptr1(enc_r[has_intersection], first_query, comprehensive_mask) # value 1 means invalid
            b = logits.size(0)
            all_neg_inf = (logits == float('-inf')).all(dim=1).to(self.device)
            idxs = torch.full((b,), self.max_degree, dtype=torch.long).to(self.device)

            log_prob = torch.zeros(b).to(self.device)
            
            valid_mask = ~all_neg_inf
            if valid_mask.any():
                valid_logits = logits[valid_mask]
                distr = Categorical(logits=valid_logits)
                if deterministic:
                    samples = valid_logits.max(dim=-1)[1]
                else:
                    samples = distr.sample()  
                idxs[valid_mask] = samples
                log_prob[valid_mask] = distr.log_prob(samples)
            # update first_idx
            intersected_indices = batch_range[has_intersection]
            new_first_idx[has_intersection] = idxs
            
            # generate new second_idx using new first_idx
            if has_intersection.any():
                # update q3 and second_query
                new_q3[intersected_indices] = encodings[intersected_indices, new_first_idx[has_intersection]]
                second_query = self.relu(first_q[has_intersection].add_(self.q_l3(new_q3[intersected_indices])))
                
                # compute the new second_idx
                unvisited = (~self.visited[has_intersection]).repeat(1, 2)
                logits = self.ptr2(enc_xy[has_intersection], second_query, unvisited)
                b = logits.size(0)
                all_neg_inf = (logits == float('-inf')).all(dim=1).to(self.device)

                log_prob = torch.zeros(b).to(self.device)

                # new logits
                valid_mask = ~all_neg_inf
                idxs_to_update = has_intersection.nonzero(as_tuple=True)[0][valid_mask]
                if valid_mask.any():
                    valid_logits = logits[valid_mask]
                    distr = Categorical(logits=valid_logits)
                    if deterministic:
                        samples = valid_logits.max(dim=-1)[1]
                    else:
                        samples = distr.sample()
                    new_idxs[idxs_to_update] = samples
                    log_prob[valid_mask] = distr.log_prob(samples)
                    # update second_idx
                    second_idx[idxs_to_update] = new_idxs[idxs_to_update] % total_length

            sec_dir[has_intersection] = torch.div(second_idx[has_intersection], self.total_length, rounding_mode='floor')
            fir_dir[has_intersection] = 1 - sec_dir[has_intersection]
            
            with torch.no_grad():
                x_idx[has_intersection] = new_first_idx[has_intersection] * sec_dir[has_intersection] + second_idx[has_intersection] * fir_dir[has_intersection]
                y_idx[has_intersection] = new_first_idx[has_intersection] * fir_dir[has_intersection] + second_idx[has_intersection] * sec_dir[has_intersection]
                first_points[has_intersection] = inputs_with_obstacle_tensor[
                                                        intersected_indices,
                                                        x_idx[has_intersection]
                                                    ]
                second_points[has_intersection] = inputs_with_obstacle_tensor[
                                                        intersected_indices,
                                                        y_idx[has_intersection]
                                                    ]
        # re-initialize the mask
        self.unactivated_obstacle = self.init_unactivated_obstacle.detach().clone()
        self.visited.scatter_(1, new_first_idx.unsqueeze(-1), True)
        
        return new_first_idx, second_idx, new_idxs, x_idx, y_idx, new_q3

class Critic(pl.LightningModule):
    def __init__(self, max_degree):
        super().__init__()
        self.save_hyperparameters()
        self.max_degree = max_degree

        # embedder args
        self.d_input = 2
        self.d_model = 128

        # encoder args
        self.num_stacks = 3
        self.num_heads = 16
        self.d_k = 16
        self.d_v = 16
        self.d_inner = 512
        self.d_unit = 256

        self.crit_embedder = Embedder(max_degree, self.d_input, self.d_model)
        # self.crit_embedder_ob_points = Embedder(max_degree, self.d_input, self.d_model)
        self.crit_encoder = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        # self.crit_encoder_ob_points = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        # self.crit_encoder_multi = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        # self.crit_embedder_ob = nn.Linear(3*self.d_model, self.d_model, bias = False)
        
        self.glimpse = Glimpse(self.d_model, self.d_unit)
        self.critic_l1 = nn.Linear(self.d_model, self.d_unit)
        self.critic_l2 = nn.Linear(self.d_unit, 1)
        # self.critic_l1_oa_length = nn.Linear(self.d_model, self.d_unit)
        # self.critic_l2_oa_length = nn.Linear(self.d_unit, 1)
        # self.critic_l1_overlap = nn.Linear(self.d_model, self.d_unit)
        # self.critic_l2_overlap = nn.Linear(self.d_unit, 1)
        self.relu = nn.ReLU()

        # self.to(device)
        self.train()

    def forward(self, inputs, input_obstacle_points_batch=None, max_obstacle=0, deterministic=False):
        if max_obstacle == 0:
            return self.forward_without_obstacle(inputs, deterministic)
        # else:
        #     return self.forward_with_obstacle(inputs, input_obstacle_points_batch, deterministic)
            
    def forward_without_obstacle(self, inputs, deterministic=False):
        inputs_tensor = torch.as_tensor(inputs).to(dtype=torch.float, device=self.device)
        # create mask for inputs_tensor
        with torch.no_grad():
            inputs_mask = (inputs_tensor[:, :, 0] == -1) # only check the first axis
        critic_encode = self.crit_encoder(self.crit_embedder(inputs_tensor)*(~inputs_mask).float().unsqueeze(-1), None)
        glimpse = self.glimpse(critic_encode)
        critic_inner = self.relu(self.critic_l1(glimpse))
        predictions = self.relu(self.critic_l2(critic_inner)).squeeze(-1)

        return predictions, None
        
    # def forward_with_obstacle(self, inputs, input_obstacle_points_batch, deterministic=False):
    #     inputs_tensor = torch.as_tensor(inputs).to(dtype=torch.float, device=self.device)
    #     # create mask for inputs_tensor
    #     with torch.no_grad():
    #         inputs_mask = (inputs_tensor[:, :, 0] == -1) # only check the first axis
    #         obstacle_mask = (input_obstacle_points_tensor[:, :, 0] == -1)
    #     input_obstacle_points_tensor = torch.as_tensor(input_obstacle_points_batch).to(dtype=torch.float, device=self.device)
    #     crit_embedings_ob_points = self.crit_embedder_ob_points(input_obstacle_points_tensor)*(~obstacle_mask).float().unsqueeze(-1)
    #     crit_embedings_ob = torch.cat((crit_embedings_ob_points[:,::2,:], crit_embedings_ob_points[:,1::2,:], crit_embedings_ob_points[:,::2,:]*crit_embedings_ob_points[:,1::2,:]), axis=2)
    #     crit_embedings_ob = self.crit_embedder_ob(crit_embedings_ob)

    #     critic_encode = self.crit_encoder(self.crit_embedder(inputs_tensor)*(~inputs_mask).float().unsqueeze(-1), None, crit_embedings_ob)
    #     critic_encode_ob_points = self.crit_encoder_ob_points(crit_embedings_ob_points, None, crit_embedings_ob)
        
    #     critic_encode_multi = self.crit_encoder_multi(torch.cat([critic_encode, critic_encode_ob_points], axis = 1), None)
    #     glimpse = self.glimpse(critic_encode_multi)
    #     # critic_inner = self.relu(self.critic_l1(glimpse))
    #     # predictions = self.relu(self.critic_l2(critic_inner)).squeeze(-1)
    #     critic_inner_oa_length = self.relu(self.critic_l1_oa_length(glimpse))
    #     predictions_oa_length = self.relu(self.critic_l2_oa_length(critic_inner_oa_length)).squeeze(-1)
    #     critic_inner_overlap = self.relu(self.critic_l1_overlap(glimpse))
    #     predictions_overlap = self.relu(self.critic_l2_overlap(critic_inner_overlap)).squeeze(-1)
        
    #     # return predictions, predictions_oa_length, predictions_overlap
    #     return predictions_oa_length, predictions_overlap