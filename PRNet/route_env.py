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
import os
from options.test_options import TestOptions
from data import create_dataset
from models import create_model
from utils import visualizer
import torch
from torchvision import transforms

from rnd import RNDModel
import torch.optim as optim
import environment_v2
import data_generator
import testNet
import random
import math
import MST

np.set_printoptions(threshold=np.inf)
rnd = RNDModel((1, 1, 84, 84), 32*32)
forward_mse = nn.MSELoss(reduction='none')
optimizer = optim.Adam(rnd.predictor.parameters(), lr=2e-5)


def compute_intrinsic_reward(rnd, next_obs):
    target_next_feature = rnd.target(next_obs)
    predict_next_feature = rnd.predictor(next_obs)

    forward_loss = forward_mse(predict_next_feature, target_next_feature).mean(-1)
    intrinsic_reward = (target_next_feature - predict_next_feature).pow(2).sum(1) / 2
    optimizer.zero_grad()
    forward_loss.backward()

    return intrinsic_reward.item()/150


grid_parameters = data_generator.dic_generate_2(
        64, 64, 5,
        5.0, 5.0,
        [1.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
        10.0,
        10.0,
        {},
)
env = environment_v2.Environment(grid_parameters)
env.load_parameters()
env.generate_capacity()


opt = TestOptions().parse()  # get test options
# hard-code some parameters for test
opt.dataroot = './test'
opt.name = 'adaptec2_obstacle_test_64'
opt.eval = 1
opt.dataset_mode = 'critical'
opt.netG = 'resnet_9blocks'

opt.num_threads = 0  # test code only supports num_threads = 0
opt.batch_size = 1  # test code only supports batch_size = 1
opt.serial_batches = True  # disable data shuffling; comment this line if results on randomly chosen images are needed.
opt.no_flip = True  # no flip; comment this line if results on flipped images are needed.
# opt.display_id = -1   # no visdom display; the test code saves the results to a HTML file.
# dataset = create_dataset(opt)  # create a dataset given opt.dataset_mode and other options
# dataset_size = len(dataset)
model = create_model(opt)  # create a model given opt.model and other options
model.setup(opt)  # regular setup: load and print networks; create schedulers


def search_ripup_and_reroute(astar, n):
    manhattan_d = 0
    wirelength = 0
    vias = 0
    while astar.paircounter < len(env.pinpairs):
        pinpair = env.pinpairs[astar.paircounter]
        # print(pinpair)
        path, seg, via = astar.search_withdemmand(pinpair[0], pinpair[1])
        astar.pathlist.append(path)
        astar.paircounter += 1
        astar.paircounter_innet += 1
        # print(path)
        manhattan_d += abs(pinpair[1][0] - pinpair[0][0])
        manhattan_d += abs(pinpair[1][1] - pinpair[0][1])
        wirelength += seg
        vias += via

    astar.passed = np.zeros_like(env.capacity)
    astar.netcouter = 0
    astar.paircounter = 0
    astar.paircounter_innet = 0
    astar.reset_pass()
    newwl = 0
    newvias = 0
    pinnum = n
    seg,via = astar.checkpath_net(pinnum)
    newwl += seg
    newvias += via
    # print(astar.newpathlist)
    return newwl, astar.newpathlist


def get_image(paths):
    inp = torch.zeros((3, 64, 64))
    for path in paths:
        for i in range(len(path)-1):
            x1, y1, z1 = path[i]
            x2, y2, z2 = path[i+1]
            if z1 == z2:
                if x1 == x2:
                    l1, l2 = min(y1, y2), max(y1, y2)
                    inp[0, x1, l1:(l2+1)] = 1
                else:
                    l1, l2 = min(x1, x2), max(x1, x2)
                    inp[0, l1:(l2+1), y1] = 1
    return transforms.Normalize((0.5,0.5,0.5),(0.5,0.5,0.5))(inp).unsqueeze_(0)


def get_training_set(nets, result):
    training_set = []
    cap_x = torch.zeros((64, 64))
    cap_y = torch.zeros((64, 64))
    for net in nets:
        # print(net)
        real_A = get_input(net, result, cap_x, cap_y)
        place = []
        for pin in net:
            x, y = result[pin]
            place.append([x, y, 0])
        twopins = []
        for i in range(len(place)):
            for j in range(i + 1, len(place)):
                twopins.append([tuple(place[i]), tuple(place[j])])
        mst = MST.generateMST(twopins)
        env.pinpairs = mst
        env.net_num = 1
        astar = environment_v2.Astar(env)
        _, path = search_ripup_and_reroute(astar, len(mst))
        real_B = get_image(path)
        training_set.append([real_A, real_B])
    return training_set


def train_model(training_set):
    index = np.random.permutation(len(training_set))
    for i in index:
        real_A, real_B = training_set[i].to('cuda:0')
        model.real_A = real_A
        model.real_B = real_B
        # model.set_input(data)
        model.optimize_parameters()


def get_input(net, place_res, h, v):
    inp = torch.zeros((3, 64, 64))
    for pin in net:
        # print(place_res[pin])
        x, y = place_res[pin]
        inp[0, x, y] = 1
    inp[1] = h
    inp[2] = v
    return transforms.Normalize((0.5,0.5,0.5),(0.5,0.5,0.5))(inp).unsqueeze_(0)


def update_cap(h, v, g):
    for i in range(len(g)):
        for j in range(len(g[0])-1):
            if g[i][j] == 1 and g[i][j+1] == 1:
                h[i][j] += 1
    for i in range(len(g[0])):
        for j in range(len(g)-1):
            if g[j][i] == 1 and g[j+1][i] == 1:
                v[j][i] += 1
    return h, v


def get_cgan_reward(nets, result, base, rate):
    wl = 0
    random.shuffle(nets)
    cap_x = torch.zeros((64, 64))
    cap_y = torch.zeros((64, 64))
    for net in nets:
        inpt = get_input(net, result, cap_x, cap_y).to('cuda:0')
        model.real_A = model.real_B = inpt
        model.forward()
        wl += model.fake_length.item()
        grid = model.fake_B[0, 0].cpu()
        cap_x, cap_y = update_cap(cap_x, cap_y, grid)
    print(wl)
    return -(wl-base) * rate, cap_x, cap_y, wl


def is_valid(x, y, size_x, size_y, ob, n):
    if -1 < x - int(size_x/2) and x + int(size_x/2) < n and -1 < y - int(size_y/2) and y + int(size_y/2) < n and int(ob[x-int(size_x/2):x+int((size_x+1)/2), y-int(size_y/2):y+int((size_y+1)/2)].sum()) < 1.0:
        return True
    return False


def search(ob, x, y, size_x, size_y, depth, n):
    if is_valid(x, y, size_x, size_y, ob, n):
        return x, y
    if depth > 9:
        return -1, -1
    elif x-1-int(size_x/2) >= 0 and is_valid(x-1, y, size_x, size_y, ob, n):
        return x-1, y
    elif x+1+int(size_x/2) < n and is_valid(x+1, y, size_x, size_y, ob, n):
        return x+1, y
    elif y-1-int(size_y/2) >= 0 and is_valid(x, y-1, size_x, size_y, ob, n):
        return x, y-1
    elif y+1+int(size_y/2) < n and is_valid(x, y+1, size_x, size_y, ob, n):
        return x, y+1
    else:
        return search(ob, x-1, y-1, size_x, size_y, depth+1, n)


def find(ob, n, size_x, size_y):
    center = [n//2, n//2]
    for i in range(n//2):
        for j in range(i):
            if is_valid(center[0]-j, center[1]-(i-j), size_x, size_y, ob, n):
                return center[0]-j, center[1]-(i-j)
            if is_valid(center[0]-j, center[1]+(i-j), size_x, size_y, ob, n):
                return center[0]-j, center[1]+(i-j)
            if is_valid(center[0]+j, center[1]-(i-j), size_x, size_y, ob, n):
                return center[0]+j, center[1]-(i-j)
            if is_valid(center[0]+j, center[1]+(i-j), size_x, size_y, ob, n):
                return center[0]+j, center[1]+(i-j)
    return -1, -1

def rule(mat):
    m, n = mat.shape
    if m == 0 or n == 0:
        return False
    else:
        return True


def edge_find(seq, ob, n, size_x, size_y):
    for node in seq:
        x = node // n
        y = node % n
        if is_valid(x, y, size_x, size_y, ob, n):
            return x, y


def edge_search(ob, x, y, size_x, size_y, n):
    if is_valid(x, y, size_x, size_y, ob, n):
        return x, y
    lx = x
    rx = n-1-x
    ly = y
    ry = n-1-y
    if lx == min(lx, rx, ly, ry):
        for i in range(lx):
            if is_valid(i, y, size_x, size_y, ob, n):
                return i, y
        if ly < ry:
            for j in range(ly):
                if is_valid(x, j, size_x, size_y, ob, n):
                    return x, j
        else:
            for j in range(ry):
                if is_valid(x, n-1-j, size_x, size_y, ob, n):
                    return x, n-1-j
    elif rx == min(lx, rx, ly, ry):
        for i in range(rx):
            if is_valid(n-1-i, y, size_x, size_y, ob, n):
                return n-1-i, y
        if ly < ry:
            for j in range(ly):
                if is_valid(x, j, size_x, size_y, ob, n):
                    return x, j
        else:
            for j in range(ry):
                if is_valid(x, n-1-j, size_x, size_y, ob, n):
                    return x, n-1-j
    elif ly == min(lx, rx, ly, ry):
        for j in range(ly):
            if is_valid(x, j, size_x, size_y, ob, n):
                return x, j
        if lx < rx:
            for i in range(lx):
                if is_valid(i, y, size_x, size_y, ob, n):
                    return i, y
        else:
            for i in range(rx):
                if is_valid(n-1-i, y, size_x, size_y, ob, n):   
                    return n-1-i, y
    elif ry == min(lx, rx, ly, ry):
        for j in range(ry):
            if is_valid(x, n-1-j, size_x, size_y, ob, n):
                return x, n-1-j
        if lx < rx:
            for i in range(lx):
                if is_valid(i, y, size_x, size_y, ob, n):
                    return i, y
        else:
            for i in range(rx):
                if is_valid(n-1-i, y, size_x, size_y, ob, n):
                    return n-1-i, y
    return -1, -1


def hpwl(r, x, n, base, rate):
    wl = 0
    con = np.zeros((n, n))
    for net in x:
        left = n-1
        right = 0
        up = n-1
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
    return -(wl-base) * rate, np.mean(con[:64])


def find_index(ava, act):
    if ava[act] == 1:
        return act
    n = len(ava)
    for i in range(1, 20):
        if act + i < n and ava[act+i] == 1:
            return act+i
        if 0 < act - i and ava[act-i] == 1:
            return act-i
    return -1


def search_index(ava, order):
    n = len(ava)
    for i in range(n):
        if ava[order[i]] == 1:
            return order[i]


def update_l(lam, n):
    return 1-math.e**(-0.006*n)


class Placememt():
    def __init__(self, name, num, base1, base2, rate1, rate2, grid_size=64):
        self.n = grid_size
        self.cell_num = num
        self.base1 = base1
        self.base2 = base2
        self.rate1 = rate1
        self.rate2 = rate2
        self.action_space = Discrete(self.n * self.n)
        self.obs_space = (1, 84, 84)
        self.obs = torch.zeros((1, 1, self.n, self.n))
        self.results = []
        self.best = 100000
        self.f = open('./result/result_' + name + '.txt', 'w')
        self.wl = 0
        self.capacity = torch.zeros((1, 1, self.n, self.n))

        f = open('./data/mix_edges_' + str(num) + '.dat', "r")
        for line in f:
            self.net = eval(line)
        f = open('./data/sizes_' + str(num) + '.dat', "r")
        for line in f:
            self.size = eval(line)
        self.init = []
        self.steps = num + len(self.net)
        self.avail = [1 for i in range(len(self.net))]
        
        def take_len(elem):
            return len(elem[1])
        self.order = [index for index, value in sorted(enumerate(self.net), key=take_len)]

        self.seed()
        self.lamb = 0
        self.seq = []
        self.place = []
        self.sequ = self.sequence()
        self.iter = 0
    
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
    
    def reset(self):
        self.obs = torch.zeros((1, 1, self.n, self.n))
        for node in self.init:
            x, y, size_x, size_y = node[0], node[1], node[2], node[3]
            self.obs[0, 0, x-int(size_x/2):x+int((size_x+1)/2), y-int(size_y/2):y+int((size_y+1)/2)] = 1
        self.capacity = torch.zeros((1, 1, self.n, self.n))
        self.seq = []
        self.place = []
        self.avail = [1 for i in range(len(self.net))]
        return self.obs

    def transform(self, x):
        up = nn.Upsample(size=84, mode='bilinear', align_corners=False)
        x = up(x)
        if torch.max(x) > 0:
            x = x / torch.max(x)
        # print(x)
        return x
    
    def sequence(self):
        n = self.n
        mat = np.array([[n*i+j for j in range(n)] for i in range(n)])
        res = list()
        while True:
            res = res + mat[0, :].tolist()
            mat = np.delete(mat, 0, 0)
            if not rule(mat):
                break
            res = res + mat[:, -1].tolist()
            mat = np.delete(mat, -1, 1)
            if not rule(mat):
                break
            res = res + list(reversed(mat[-1, :].tolist()))
            mat = np.delete(mat, -1, 0)
            if not rule(mat):
                break
            res = res + list(reversed(mat[:, 0].tolist()))
            mat = np.delete(mat, 0, 1)
            if not rule(mat):
                break
        return res
    
    def step(self, action, index):
        if index < self.cell_num:
            size_x = self.size[index][0]
            size_y = self.size[index][1]
            x = action // self.n
            y = action % self.n
            x, y = search(self.obs[0, 0], x, y, size_x, size_y, 0, self.n)
            if x == -1 or y == -1:
                x, y = find(self.obs[0, 0], self.n, size_x, size_y)
            self.obs[0, 0, x-int(size_x/2):x+int((size_x+1)/2), y-int(size_y/2):y+int((size_y+1)/2)] = 1
            self.results.append([int(x), int(y)])
            if len(self.results) == self.cell_num:
                net_feature = torch.zeros(len(self.net), 11)
                for i in range(len(self.net)):
                    n_i = self.net[i]
                    for j in range(len(n_i)):
                        x, y = self.results[n_i[j]]
                        net_feature[i, 2*j] = x
                        net_feature[i, 2*j+1] = y
                torch.save(net_feature, './data/net_feature_' + str(self.cell_num) + '.pt')
            self.place = self.results.copy()
            obs = self.transform(self.obs)
            done = False
            reward = compute_intrinsic_reward(rnd, obs / 255.0)
            return obs, done, torch.FloatTensor([[reward]])
        
        else:
            i_net = find_index(self.avail, action)
            if i_net == -1:
                i_net = search_index(self.avail, self.order)
            net_feature = torch.load('./data/net_feature_' + str(self.cell_num) + '.pt')
            net_feature[i_net, 10] = 1
            torch.save(net_feature, './data/net_feature_' + str(self.cell_num) + '.pt')
            self.avail[i_net] = 0
            net_i = self.net[i_net]
            self.seq.append(net_i)

            for pin in net_i:
                x, y = self.results[pin]
                self.capacity[0, 0, x, y] += 1

            cap = self.transform(self.capacity)
            self.results.append([int(action)])

            if len(self.results) == self.steps:
                done = True
                self.iter += 1
                if self.iter % 20 == 0:
                    training_set = get_training_set(self.net, self.results)
                    train_model(training_set)
                wl1, cp_x, cp_y, wl = get_cgan_reward(self.net, self.results, self.base1, self.rate1)
                wl2, con = hpwl(self.results, self.net, self.n, self.base2, self.rate2)
                # print(wl1, wl2, self.lamb)
                reward = self.lamb * wl1 + (1 - self.lamb) * wl2
                self.lamb = update_l(self.lamb, self.iter)
                if wl < self.best:
                    self.best = wl
                    self.f.write(str(self.obs))
                    self.f.write(str(self.results))
                    self.f.write('\n')
                    self.f.write(str(wl))
                    self.f.write('\n')
                    self.f.write(str(con))
                    self.f.write('\n')
                    torch.save([cp_x, cp_y], './data/cap_' + str(self.steps) + '.pt')
                self.results = []
            else:
                done = False
                reward = compute_intrinsic_reward(rnd, cap)
            return cap, done, torch.FloatTensor([[reward]])


def route_envs(name, num, base1, base2, rate1, rate2):
    return Placememt(name, num, base1, base2, rate1, rate2)
