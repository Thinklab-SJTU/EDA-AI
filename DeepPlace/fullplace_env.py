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
import dreamplace.configure as configure
import Params
import PlaceDB
import NonLinearPlace
import pdb

from rnd import RNDModel
import torch.optim as optim


np.set_printoptions(threshold=np.inf)
rnd = RNDModel((1, 1, 84, 84), 32*32)
forward_mse = nn.MSELoss(reduction='none')
optimizer = optim.Adam(rnd.predictor.parameters(), lr=2e-6)

def compute_intrinsic_reward(rnd, next_obs):
    target_next_feature = rnd.target(next_obs)
    predict_next_feature = rnd.predictor(next_obs)

    forward_loss = forward_mse(predict_next_feature, target_next_feature).mean(-1)
    intrinsic_reward = (target_next_feature - predict_next_feature).pow(2).sum(1) / 2
    optimizer.zero_grad()
    forward_loss.backward()

    return intrinsic_reward.item()/200


def place(params):
    """
    @brief Top API to run the entire placement flow. 
    @param params parameters 
    """

    assert (not params.gpu) or configure.compile_configurations["CUDA_FOUND"] == 'TRUE', \
            "CANNOT enable GPU without CUDA compiled"

    np.random.seed(params.random_seed)
    # read database
    placedb = PlaceDB.PlaceDB()
    placedb(params)

    # solve placement
    placer = NonLinearPlace.NonLinearPlace(params, placedb)
    metrics = placer(params, placedb)
    result = metrics[-3][0]

    # write placement solution
    path = "%s/%s" % (params.result_dir, params.design_name())
    if not os.path.exists(path):
        os.system("mkdir -p %s" % (path))
    gp_out_file = os.path.join(
        path,
        "%s.gp.%s" % (params.design_name(), params.solution_file_suffix()))
    placedb.write(params, gp_out_file)

    # call external detailed placement
    # TODO: support more external placers, currently only support
    # 1. NTUplace3/NTUplace4h with Bookshelf format
    # 2. NTUplace_4dr with LEF/DEF format
    if params.detailed_place_engine and os.path.exists(
            params.detailed_place_engine):
        logging.info("Use external detailed placement engine %s" %
                     (params.detailed_place_engine))
        if params.solution_file_suffix() == "pl" and any(
                dp_engine in params.detailed_place_engine
                for dp_engine in ['ntuplace3', 'ntuplace4h']):
            dp_out_file = gp_out_file.replace(".gp.pl", "")
            # add target density constraint if provided
            target_density_cmd = ""
            if params.target_density < 1.0 and not params.routability_opt_flag:
                target_density_cmd = " -util %f" % (params.target_density)
            cmd = "%s -aux %s -loadpl %s %s -out %s -noglobal %s" % (
                params.detailed_place_engine, params.aux_input, gp_out_file,
                target_density_cmd, dp_out_file, params.detailed_place_command)
            logging.info("%s" % (cmd))
            # tt = time.time()
            os.system(cmd)
            # logging.info("External detailed placement takes %.2f seconds" %
            #              (time.time() - tt))

            if params.plot_flag:
                # read solution and evaluate
                placedb.read_pl(params, dp_out_file + ".ntup.pl")
                iteration = len(metrics)
                pos = placer.init_pos
                pos[0:placedb.num_physical_nodes] = placedb.node_x
                pos[placedb.num_nodes:placedb.num_nodes +
                    placedb.num_physical_nodes] = placedb.node_y
                hpwl, density_overflow, max_density = placer.validate(
                    placedb, pos, iteration)
                logging.info(
                    "iteration %4d, HPWL %.3E, overflow %.3E, max density %.3E"
                    % (iteration, hpwl, density_overflow, max_density))
                placer.plot(params, placedb, iteration, pos)
        elif 'ntuplace_4dr' in params.detailed_place_engine:
            dp_out_file = gp_out_file.replace(".gp.def", "")
            cmd = "%s" % (params.detailed_place_engine)
            for lef in params.lef_input:
                if "tech.lef" in lef:
                    cmd += " -tech_lef %s" % (lef)
                else:
                    cmd += " -cell_lef %s" % (lef)
            cmd += " -floorplan_def %s" % (gp_out_file)
            cmd += " -verilog %s" % (params.verilog_input)
            cmd += " -out ntuplace_4dr_out"
            cmd += " -placement_constraints %s/placement.constraints" % (
                os.path.dirname(params.verilog_input))
            cmd += " -noglobal %s ; " % (params.detailed_place_command)
            cmd += "mv ntuplace_4dr_out.fence.plt %s.fense.plt ; " % (
                dp_out_file)
            cmd += "mv ntuplace_4dr_out.init.plt %s.init.plt ; " % (
                dp_out_file)
            cmd += "mv ntuplace_4dr_out %s.ntup.def ; " % (dp_out_file)
            cmd += "mv ntuplace_4dr_out.ntup.overflow.plt %s.ntup.overflow.plt ; " % (
                dp_out_file)
            cmd += "mv ntuplace_4dr_out.ntup.plt %s.ntup.plt ; " % (
                dp_out_file)
            if os.path.exists("%s/dat" % (os.path.dirname(dp_out_file))):
                cmd += "rm -r %s/dat ; " % (os.path.dirname(dp_out_file))
            cmd += "mv dat %s/ ; " % (os.path.dirname(dp_out_file))
            # logging.info("%s" % (cmd))
            # tt = time.time()
            os.system(cmd)
            # logging.info("External detailed placement takes %.2f seconds" %
            #              (time.time() - tt))
        else:
            logging.warning(
                "External detailed placement only supports NTUplace3/NTUplace4dr API"
            )
    elif params.detailed_place_engine:
        logging.warning(
            "External detailed placement engine %s or aux file NOT found" %
            (params.detailed_place_engine))

    return result


def write(res):
    dic = np.load('./DeepPlace/data/3_dic.npy', allow_pickle=True).item()
    f = open("./benchmarks/ispd2005/adaptec3/adaptec3.pl", "w")
    with open("./DeepPlace/data/adaptec3.pl", "r") as f2:
        for line in f2:
            line = line.strip()
            l = line.split()
            if line and l[0][0] == 'o':
                num = int(l[0].lstrip('o'))
                if num - 450927 in dic.keys():
                    index = dic[num - 450927]
                    pos = res[index]
                    x = int(pos[0] / 32 * 22653)
                    y = int(pos[1] / 32 * 23122)
                    l[1] = str(x)
                    l[2] = str(y)
                    line = '\t'.join(l)

            f.write(line)
            f.write('\n')


def new_cal_re(res, params):
    write(res)
    r = place(params)
    wl = float(r[0].hpwl.data)
    overf = float(r[0].overflow.data)
    reward = -2 * (wl - 2.4e8) * 1e-6 - overf * 20
    print(reward, wl)
    return reward


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


def is_valid(x, y):
    if -1 < x < 32 and -1 < y < 32:
        return True
    return False


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


class Placememt():
    def __init__(self, grid_size=32, num_cell=710):
        self.n = grid_size
        self.steps = num_cell
        self.action_space = Discrete(self.n * self.n)
        self.obs_space = (1, 84, 84)
        self.obs = torch.zeros((1, 1, self.n, self.n))
        self.results = []
        self.best = -500
        self.f = open("./DeepPlace/result/result.txt", 'w')

        f = open("./DeepPlace/data/n_edges_710.dat", "r")
        for line in f:
            self.net = eval(line)
        self.seed()
        logging.root.name = 'DREAMPlace'
        self.params = Params.Params()

        # load parameters
        add = "test/ispd2005/adaptec3.json"
        self.params.load(add)
        os.environ["OMP_NUM_THREADS"] = "%d" % (self.params.num_threads)
    
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
            reward = new_cal_re(self.results, self.params)
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


def fullplace_envs():
    return Placememt()
