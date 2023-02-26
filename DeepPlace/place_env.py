from gym.spaces import Discrete
import torch
import torch.nn as nn
import numpy as np
from gym.utils import seeding
import os
import sys
import logging
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if root_dir not in sys.path:
    sys.path.append(root_dir)

from rnd import RNDModel
import torch.optim as optim


np.set_printoptions(threshold=np.inf)
rnd = RNDModel((1, 1, 84, 84), 32*32)
forward_mse = nn.MSELoss(reduction='none')
optimizer = optim.Adam(rnd.predictor.parameters(), lr=5e-6)


def compute_intrinsic_reward(rnd, next_obs):
    target_next_feature = rnd.target(next_obs)
    predict_next_feature = rnd.predictor(next_obs)

    forward_loss = forward_mse(predict_next_feature, target_next_feature).mean(-1)
    intrinsic_reward = (target_next_feature - predict_next_feature).pow(2).sum(1) / 2
    optimizer.zero_grad()
    forward_loss.backward()

    return intrinsic_reward.item()/100


def is_valid(x, y):
    if -1 < x < 32 and -1 < y < 32:
        return True
    return False


def search(ob, x, y, depth, n):
    if ob[x, y] < 1.0:
        return x, y
    if depth > 7:
        return -1, -1
    elif x-1 >= 0 and ob[x-1, y] < 1.0:
        return x-1, y
    elif x+1 < n and ob[x+1, y] < 1.0:
        return x+1, y
    elif y-1 >= 0 and ob[x, y-1] < 1.0:
        return x, y-1
    elif y+1 < n and ob[x, y+1] < 1.0:
        return x, y+1
    else:
        return search(ob, x-1, y-1, depth+1, n)


def find(ob, n):
    center = [n//2, n//2]
    for i in range(n):
        for j in range(i):
            if is_valid(center[0]-j, center[1]-(i-j)) and ob[center[0]-j, center[1]-(i-j)] < 1.0:
                return center[0]-j, center[1]-(i-j)
            if is_valid(center[0]-j, center[1]+(i-j)) and ob[center[0]-j, center[1]+(i-j)] < 1.0:
                return center[0]-j, center[1]+(i-j)
            if is_valid(center[0]+j, center[1]-(i-j)) and ob[center[0]+j, center[1]-(i-j)] < 1.0:
                return center[0]+j, center[1]-(i-j)
            if is_valid(center[0]+j, center[1]+(i-j)) and ob[center[0]+j, center[1]+(i-j)] < 1.0:
                return center[0]+j, center[1]+(i-j)


def cal_re(r, x):
    wl = 0
    con = np.zeros((32, 32))
    for net in x:
        left = 31
        right = 0
        up = 31
        down = 0
        for i in net:
            left = min(left, r[i][1])
            right = max(right, r[i][1])
            up = min(up, r[i][0])
            down = max(down, r[i][0])
        wn = int(right-left+1)
        hn = int(down-up+1)
        dn = (wn+hn) / (wn*hn)
        con[up:down+1, left:right+1] += dn
        wl += wn + hn
    con = list(con.flatten())
    con.sort(reverse=True)
    return (-np.mean(con[:32]) - (wl-34000)*0.1)*0.2

class Placememt():
    def __init__(self, grid_size=32, num_cell=710):
        self.n = grid_size
        self.steps = num_cell
        self.action_space = Discrete(self.n * self.n)
        self.obs_space = (1, 84, 84)
        self.obs = torch.zeros((1, 1, self.n, self.n))
        self.results = []
        self.best = -500
        self.f = open("./result/result.txt", 'w')

        f = open("./data/n_edges_710.dat", "r")
        for line in f:
            self.net = eval(line)
        self.seed()

    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    def reset(self):
        self.obs = torch.zeros((1, 1, self.n, self.n))
        return self.obs
    
    def to(self, device):
        self.obs = self.obs.to(device)
        
    def transform(self, x):
        up = nn.Upsample(size=84, mode='bilinear', align_corners=False)
        return up(x)*255
    
    def step(self, action):
        x = action // self.n
        y = action % self.n
        x, y = search(self.obs[0, 0], x, y, 0, self.n)
        if x == -1 or y == -1:
            x, y = find(self.obs[0, 0], self.n)
        self.obs[0, 0, x, y] = 1
        self.results.append([int(x), int(y)])
        obs = self.transform(self.obs)

        if len(self.results) == self.steps:
            done = True
            reward = cal_re(self.results, self.net)
            if reward > self.best:
                self.best = reward
                self.f.write(str(self.obs))
                self.f.write(str(self.results))
                self.f.write('\n')
                self.f.write(str(reward))
                self.f.write('\n')
            self.results = []
        else:
            done = False
            reward = compute_intrinsic_reward(rnd, obs / 255.0)
        return obs, done, torch.FloatTensor([[reward]])


def place_envs():
    return Placememt()
