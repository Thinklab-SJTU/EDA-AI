import numpy as np
from utils import *
from queue import PriorityQueue
import data_generator
import random
from heapq import *

MAX_STEP = 120
MOVING_RATE = 0.5
MAX_EPISODE=500
TH=20

def manhattan_distance(pos1, pos2):
    res = 0
    for i in range(3):
        res += abs(pos1[i] - pos2[i])
    return res


class Environment:
    def __init__(self, grid_parameters):
        self.gridParameters = grid_parameters
        self.width = 0
        self.height = 0
        self.layers = 0
        self.vertical_capacity = None
        self.horizontal_capacity = None
        self.layer_capacity = None
        self.capacity = None
        self.capacity_cache = None
        self.capacity_prev= None

        self.specified_capacity = None  # ([x,y,z],[c1,c2,c3,c4,c5,c6])

        self.current_state = None
        self.target_state = None

        self.current_route = []
        self.best_route=[]

        self.current_reward = 0
        self.best_reward=0
        self.reward_log = []

        self.original_capacity = None
        self.steplength = None
        self.passed = None  # 记录已经布过的线，在同一个Net中重复经过的时候不需要更新capacity
        self.passed_cache=None
        self.step_counter = 0

        self.pinpairs = None
        self.pinpairs_innet=None
        self.pair_counter = 0
        self.pair_counter_innet = 0  # 记录当前net做完了多少个pair
        self.net_num = None  # 记录每个Net有多少个pair
        self.net_counter = 0
        self.pair_num = 0

        self.episode_counter = 0

        self.route = []
        self.route_log = []

        self.route_window = []
        self.padding = 2

        self.learning_counter_inpair=0
        self.cvg_counter=0
        self.converged=False
        self.failed_pairs=[]
        self.failed_net=[]
        return

    def load_parameters(self):
        gp = self.gridParameters
        grid_size = gp['gridSize']
        self.width = grid_size[0]
        self.height = grid_size[1]
        self.layers = grid_size[2]
        self.vertical_capacity = []
        for i in range(self.layers):
            self.vertical_capacity.append(int(gp['verticalCapacity'][i]))
        self.horizontal_capacity = []
        for i in range(self.layers):
            self.horizontal_capacity.append(int(gp['horizontalCapacity'][i]))
        self.layer_capacity = []
        for i in range(self.layers):
            self.layer_capacity.append(20)
        self.steplength = [int(gp['tileWidth']), int(gp['tileHeight'])]
        if 'netInfo' in gp:
            self.pinpairs = []
            pairs_by_net = generatePairs(gp)
            self.net_num = []
            for net in pairs_by_net:
                self.net_num.append(len(net))
                for pair in net:
                    self.pinpairs.append(pair)

            self.pair_num = len(self.pinpairs)

    def state_to_input(self, state):
        # print(state)
        x, y, z = state[:3]
        state = np.array(state[:3])
        capacity = np.squeeze(self.capacity[x, y, z, :])
        distance = np.array(self.target_state[:3]) - state
        net_input = np.concatenate((state, capacity, distance), axis=0).reshape(1, -1)
        return net_input

    def generate_capacity(self):
        self.capacity = np.zeros((self.width, self.height, self.layers, 6))
        # 实际上三个len都是一样的 都是层数
        for i in range(self.layers):
            self.capacity[:, :, i, 0] = self.capacity[:, :, i, 1] = self.horizontal_capacity[i]

        for i in range(self.layers):
            self.capacity[:, :, i, 2] = self.capacity[:, :, i, 3] = self.vertical_capacity[i]

        for i in range(self.layers):
            self.capacity[:, :, i, 4] = self.capacity[:, :, i, 5] = self.layer_capacity[i]

        if self.specified_capacity:
            for i in range(len(self.specified_capacity)):
                x, y, z = self.specified_capacity[i][0]
                self.capacity[x, y, z, :] = self.specified_capacity[i][1]

        self.capacity[:, :, -1, 4] = 0
        self.capacity[:, :, 0, 5] = 0

        self.capacity[:, 0, :, 3] = 0
        self.capacity[:, self.height - 1, :, 2] = 0

        self.capacity[0, :, :, 1] = 0
        self.capacity[self.width - 1, :, :, 0] = 0

        self.original_capacity = self.capacity
        self.passed = np.zeros_like(self.capacity)
        self.passed_cache=self.passed.copy()
        self.passed_prev=self.passed.copy()

        self.capacity_cache = self.capacity.copy()
        self.capacity_prev=self.capacity.copy()

        return

    def step(self, action):
        state = self.current_state
        x, y, z = state[0], state[1], state[2]
        #real_x, real_y = state[3], state[4]
        moved = True
        reward = 0
        next_state = state.copy()
        if action == 0 and self.capacity_cache[x, y, z, 0] > 0 and x + 1 <= self.route_window[1]:
            next_state[0] += 1
            #next_state[3] += self.steplength[0]

        elif action == 1 and self.capacity_cache[x, y, z, 1] > 0 and x - 1 >= self.route_window[0]:
            next_state[0] -= 1
            #next_state[3] -= self.steplength[0]

        elif action == 2 and self.capacity_cache[x, y, z, 2] > 0 and y + 1 <= self.route_window[3]:
            next_state[1] += 1
            #next_state[4] += self.steplength[1]

        elif action == 3 and self.capacity_cache[x, y, z, 3] > 0 and y - 1 >= self.route_window[2]:
            next_state[1] -= 1
            #next_state[4] -= self.steplength[1]

        elif action == 4 and self.capacity_cache[x, y, z, 4] > 0:
            next_state[2] += 1

        elif action == 5 and self.capacity_cache[x, y, z, 5] > 0:
            next_state[2] -= 1
        else:
            moved = False

        if moved:
            if self.passed_cache[x, y, z, action] == 0:
                self.passed_cache[x, y, z, action] = 1

                # update capacity
                self.update_capacity_cache(state, action, next_state)
            self.route.append(self.current_state.copy())

        self.current_state = next_state
        self.step_counter += 1

        done = False
        exceeded = False

        if moved:
            reward = -MOVING_RATE
        else:
            reward = -1

        if self.current_state[:3] == self.target_state[:3]:
            #print("#", self.current_state, self.target_state)
            done = True
            reward = MAX_STEP

        elif self.step_counter == MAX_STEP:
            done = True
            exceeded = True

        self.current_reward += reward

        if done:
            if not exceeded:
                if self.current_reward>self.best_reward:
                    self.capacity=self.capacity_cache.copy()
                    self.passed=self.passed_cache.copy()
                    self.best_route=self.route.copy()
                    self.best_reward=self.current_reward
                self.cvg_counter+=1
                self.learning_counter_inpair+=1
                if self.cvg_counter==TH:
                    self.converged=True

            else:
                self.cvg_counter=0
                self.learning_counter_inpair+=1



        return next_state, reward, done, [exceeded]


    def output(self):
        rsum = 0
        succeeded = 0
        for i in range(self.pair_num):
            route = self.route_log[i]
            reward = self.reward_log[i]
            if reward >= 0:
                rsum += reward
                succeeded += 1
            print(reward, route)
        print("total reward ", rsum, "routed pairs ", succeeded)

    def reset(self):

        if self.learning_counter_inpair==MAX_EPISODE or self.converged==True:
            if self.converged:
                self.capacity_prev = self.capacity.copy()
                self.passed_prev=self.passed.copy()

            else:
                self.failed_pairs.append(self.pinpairs[self.pair_counter])

            self.pair_counter += 1
            self.pair_counter_innet += 1

            if self.pair_counter_innet == self.net_num[self.net_counter]:  # 一个Net中的Pinpair全部做完
                self.net_counter += 1
                self.pair_counter_innet = 0
                self.passed = np.zeros_like(self.capacity)
                self.passed_prev=self.passed.copy()
                self.passed_cache=self.passed.copy()

            if self.converged:

                self.route_log.append(self.best_route)
                self.best_reward = round(self.best_reward)
                self.reward_log.append(self.best_reward)
            else:
                self.route_log.append([])
                self.reward_log.append(0)

            self.best_reward=0
            self.learning_counter_inpair = 0
            self.cvg_counter=0
            self.converged=False

            if self.pair_counter == self.pair_num:
                return [],True

        self.capacity_cache = self.capacity_prev.copy()
        self.passed_cache=self.passed_prev.copy()

        self.current_state = list(self.pinpairs[self.pair_counter][0])
        self.target_state = list(self.pinpairs[self.pair_counter][1])

        if self.learning_counter_inpair==0:
            self.generate_window(self.current_state, self.target_state)


        # print(self.current_state,self.target_state,self.route_window)

        self.route = []
        self.current_reward = 0
        self.step_counter = 0

        return self.current_state,False

    def generate_window(self, start, end):
        x1 = max(0, start[0] - self.padding)
        x2 = max(0, end[0] - self.padding)
        x3 = min(self.width - 1, start[0] + self.padding)
        x4 = min(self.width - 1, end[0] + self.padding)

        left = min(x1, x2)
        right = max(x3, x4)

        y1 = max(0, start[1] - self.padding)
        y2 = max(0, end[1] - self.padding)
        y3 = min(self.height - 1, start[1] + self.padding)
        y4 = min(self.height - 1, end[1] + self.padding)

        low = min(y1, y2)
        high = max(y3, y4)

        self.route_window = [left, right, low, high]

    def update_capacity(self, state, action, next_state):
        x, y, z = state[:3]
        self.capacity[x, y, z, action] -= 1
        reverse_action = action
        if action % 2 == 0:
            reverse_action += 1
        else:
            reverse_action -= 1
        self.capacity[next_state[0], next_state[1], next_state[2], reverse_action] -= 1

    def update_capacity_cache(self, state, action, next_state):
        x, y, z = state[:3]
        self.capacity_cache[x, y, z, action] -= 1
        reverse_action = action
        if action % 2 == 0:
            reverse_action += 1
        else:
            reverse_action -= 1
        self.capacity_cache[next_state[0], next_state[1], next_state[2], reverse_action] -= 1

    def output_capacity(self,cap):
        for i in range(self.layers):
            for j in range(self.width):
                for k in range(self.height):
                    print("%.2f"%np.mean(cap[j,k,i]),end=' ')
                print(' ')
            print(' ')
        print(' ')


class Astar:
    def __init__(self, env):
        self.width = env.width
        self.height = env.height
        self.layers = env.layers
        self.capacity = env.capacity.copy()
        self.pinpairs = env.pinpairs
        self.originalcap = env.capacity.copy()
        self.passed=np.zeros_like(self.capacity)
        self.paircounter=0
        self.netnum=env.net_num
        self.netcouter=0
        self.paircounter_innet=0
        self.demand=np.zeros_like(env.capacity)
        self.pathlist=[]
        self.newpathlist=[]


    def ripup_path(self,path):
        for i in range(len(path)-1):
            cur = path[i]
            next = path[i + 1]
            action = self.get_action_list(cur, next)
            self.passed[cur[0],cur[1],cur[2],action]-=1
            if self.passed[cur[0],cur[1],cur[2],action]==0:
                self.demand[cur[0],cur[1],cur[2],action]-=1
            reverse_action = action
            if action % 2 == 0:
                reverse_action += 1
            else:
                reverse_action -= 1
            self.passed[next[0],next[1],next[2],reverse_action]-=1
            if self.passed[next[0],next[1],next[2],reverse_action]==0:
                self.demand[next[0],next[1],next[2],reverse_action]-=1

    def checkpath_net(self,n):
        for j in range(n):
            path=self.pathlist[self.paircounter+j]
            self.update_pass_rr(path)
        total_seg=0
        total_via=0
        for j in range(n):
            path = self.pathlist[self.paircounter + j]
            seg,via=self.output_path(path)
            newpath=path
            for i in range(len(path) - 1):
                cur = path[i]
                x, y, z = cur
                next = path[i + 1]
                action = self.get_action_list(cur, next)
                if self.demand[x, y, z][action] > self.capacity[x, y, z][action]:
                    self.ripup_path(path)
                    start = path[0]
                    end = path[-1]
                    newpath, seg, via = self.search_rr(start, end)
                    break
            # print(path)
            # print(newpath)
            self.newpathlist.append(newpath)
            total_seg+=seg
            total_via+=via
       
        return total_seg,total_via

    def checkpath(self,path):
        newpath=path
        seg,via=self.output_path(path)
        updated=False
        for i in range(len(path)-1):
            cur=path[i]
            x,y,z=cur
            next=path[i+1]
            action=self.get_action_list(cur,next)
            if self.demand[x,y,z][action]>self.capacity[x,y,z][action]:
                self.ripup_path(path)
                start=path[0]
                end=path[-1]
                newpath,seg,via=self.search_rr(start,end)
                updated=True
            else:
                self.update_pass_rr(path)

        return newpath,seg,via,updated

    def reset_cap(self):
        self.capacity = self.originalcap.copy()

    def reset_pass(self):
        self.passed=np.zeros_like(self.capacity)


    def output_path(self,path):
        seg=0
        via=0
        for i in range(len(path)-1):
            cur=path[i]
            next=path[i+1]
            action=self.get_action_list(cur,next)
            if action<=3:
                seg+=1
            else:
                via+=1

        return seg,via

    def update_pass_rr(self,path):
        for i in range(len(path) - 1):
            cur = path[i]
            next = path[i + 1]
            action = self.get_action_list(cur, next)
            reverse_action = action
            if action % 2 == 0:
                reverse_action += 1
            else:
                reverse_action -= 1
            self.passed[cur[0], cur[1], cur[2], action] += 1
            self.passed[next[0], next[1], next[2], reverse_action] += 1
        return

    def update_capacity(self,cap,state,action,next_state,dif):
        x, y, z = state[:3]
        cap[x, y, z, action] +=dif
        reverse_action = action
        if action % 2 == 0:
            reverse_action += 1
        else:
            reverse_action -= 1
        cap[next_state[0], next_state[1], next_state[2], reverse_action] +=dif



    def update_path(self,info,start,end,cap,dif):
        cur = list(end)
        path = []
        path.insert(0, cur)
        seg = 0
        via = 0
        while True:
            last = path[0]

            x, y, z = info[(cur[0], cur[1], cur[2])][4:7]
            cur = [int(x), int(y), int(z)]
            action = self.get_action_list(cur, last)
            if action <= 3:
                seg += 1
            else:
                via += 1

            reverse_action = action
            if action % 2 == 0:
                reverse_action += 1
            else:
                reverse_action -= 1

            if self.passed[x, y, z][action] == 0:

                self.update_capacity(cap,cur, action, last,dif)

            self.passed[x, y, z][action] +=1
            self.passed[last[0], last[1], last[2]][reverse_action] +=1
            path.insert(0, cur)

            if cur == list(start[:3]):
                break
        return path, seg, via


    def bury_experience(self, env, network, path):
        l = len(path)
        for i in range(l - 1):
            state = env.state_to_input(path[i])
            # print(state)
            next_state = env.state_to_input(path[i + 1])
            action = self.get_action_input(state, next_state)
            reward = MAX_STEP - (l - i) * MOVING_RATE
            t = network.store_transition_astar(state, action, reward, next_state, reward)
            # print(t)

    def get_action_list(self, state, next_state):
        dif = [state[i] - next_state[i] for i in range(len(state))]
        return self.get_action(dif)

    def get_action(self,dif):
        if dif[0] == 1:
            return 0
        elif dif[0] == -1:
            return 1
        elif dif[1] == 1:
            return 2
        elif dif[1] == -1:
            return 3
        elif dif[2] == 1:
            return 4
        else:
            return 5

    def get_action_input(self, state, next_state):
        dif = (next_state - state)[0]
        return self.get_action(dif)


    def search_withcap(self, start, end,cap):
        info = {}
        openlist = []

        f = manhattan_distance(start, end)
        heappush(openlist, [f, start])

        info[tuple(start)] = [0, 1, f, 0, -1, -1, -1]  # f,g,parent

        while len(openlist) != 0:
            cur = heappop(openlist)
            x, y, z = cur[1][:3]
            info[(x, y, z)][0] = 1
            adj = [(x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)]
            for i in range(6):
                # print(adj[i])
                newx, newy, newz = adj[i]
                if newx < 0 or newx >= self.width or newy < 0 or newy >= self.height or newz < 0 or newz >= self.layers:
                    continue
                if adj[i] in info:
                    inclose, inopen, fvalue = info[(newx, newy, newz)][:3]
                    if inclose:
                        continue

                    else:  # inopen
                        if cap[x, y, z, i] > 0:
                            prevf = info[(newx, newy, newz)][2]
                            newf = info[(x, y, z)][3] + 1 + manhattan_distance(cur[1], end)
                            if newf < prevf:
                                info[(newx, newy, newz)][2] = newf
                                info[(newx, newy, newz)][3] = info[(x, y, z)][3] + 1
                                info[(newx, newy, newz)][4:7] = x, y, z
                                for item in openlist:
                                    if item[1] == [newx, newy, newz]:
                                        item[0] = newf
                                        # print(cur[1],adj[i],newf,'justified')

                else:
                    if cap[x, y, z, i] > 0:

                        g = info[(x, y, z)][3] + 1
                        h = manhattan_distance(adj[i], end)
                        info[(newx, newy, newz)] = [0, 1, g + h, g, x, y, z]

                        heappush(openlist, [g + h, adj[i]])
                        # print(cur[1], adj[i],g+h,'inserted')
                    else:
                        continue
            # print('')
            heapify(openlist)
            if (end[0], end[1], end[2]) in info:
                if info[(end[0], end[1], end[2])][0]:
                    break

        path = []
        seg = 0
        via = 0
        if (end[0], end[1], end[2]) in info:
            if info[(end[0], end[1], end[2])][0]:
                path, seg, via = self.update_path(info, start, end,cap,-1)
        return path, seg, via


    def search_withdemmand(self, start, end):
        info = {}
        openlist = []

        f = manhattan_distance(start, end)
        heappush(openlist, [f, start])

        info[tuple(start)] = [0, 1, f, 0, -1, -1, -1]  # f,g,parent

        while len(openlist) != 0:
            cur = heappop(openlist)
            x, y, z = cur[1][:3]
            info[(x, y, z)][0] = 1
            adj = [(x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)]
            for i in range(6):
                # print(adj[i])
                newx, newy, newz = adj[i]
                if newx < 0 or newx >= self.width or newy < 0 or newy >= self.height or newz < 0 or newz >= self.layers:
                    continue
                if adj[i] in info:
                    inclose, inopen, fvalue = info[(newx, newy, newz)][:3]
                    if inclose:
                        continue

                    else:  # inopen
                        if self.capacity[x,y,z,i]>0:
                            prevf = info[(newx, newy, newz)][2]
                            newf = info[(x, y, z)][3] + 1 + manhattan_distance(cur[1], end)
                            if newf < prevf:
                                info[(newx, newy, newz)][2] = newf
                                info[(newx, newy, newz)][3] = info[(x, y, z)][3] + 1
                                info[(newx, newy, newz)][4:7] = x, y, z
                                for item in openlist:
                                    if item[1] == [newx, newy, newz]:
                                        item[0] = newf
                                        # print(cur[1],adj[i],newf,'justified')

                else:
                    if self.capacity[x,y,z,i]>0:
                        g = info[(x, y, z)][3] + 1
                        h = manhattan_distance(adj[i], end)
                        info[(newx, newy, newz)] = [0, 1, g + h, g, x, y, z]

                        heappush(openlist, [g + h, adj[i]])

            # print('')
            heapify(openlist)
            if (end[0], end[1], end[2]) in info:
                if info[(end[0], end[1], end[2])][0]:
                    break

        path = []
        seg = 0
        via = 0
        if (end[0], end[1], end[2]) in info:

            if info[(end[0], end[1], end[2])][0]:
                path, seg, via = self.update_path(info, start, end,self.demand,1)
        return path, seg, via

    def search_rr(self, start, end):
        info = {}
        openlist = []

        f = manhattan_distance(start, end)
        heappush(openlist, [f, start])

        info[tuple(start)] = [0, 1, f, 0, -1, -1, -1]  # f,g,parent

        while len(openlist) != 0:
            cur = heappop(openlist)

            x, y, z = cur[1][:3]
            info[(x, y, z)][0] = 1
            adj = [(x + 1, y, z), (x - 1, y, z), (x, y + 1, z), (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)]
            for i in range(6):
                # print(adj[i])
                newx, newy, newz = adj[i]
                #print(cur,adj[i],self.demand[x,y,z,i])
                if newx < 0 or newx >= self.width or newy < 0 or newy >= self.height or newz < 0 or newz >= self.layers:
                    continue
                if adj[i] in info:
                    inclose, inopen, fvalue = info[(newx, newy, newz)][:3]
                    if inclose:
                        continue

                    else:  # inopen
                        if self.demand[x,y,z,i]<self.capacity[x,y,z,i]:
                            prevf = info[(newx, newy, newz)][2]
                            newf = info[(x, y, z)][3] + 1 + manhattan_distance(cur[1], end)
                            if newf < prevf:
                                info[(newx, newy, newz)][2] = newf
                                info[(newx, newy, newz)][3] = info[(x, y, z)][3] + 1
                                info[(newx, newy, newz)][4:7] = x, y, z
                                for item in openlist:
                                    if item[1] == [newx, newy, newz]:
                                        item[0] = newf
                                        # print(cur[1],adj[i],newf,'justified')

                else:
                    if self.demand[x, y, z, i] < self.capacity[x, y, z, i]:
                        g = info[(x, y, z)][3] + 1
                        h = manhattan_distance(adj[i], end)
                        info[(newx, newy, newz)] = [0, 1, g + h, g, x, y, z]

                        heappush(openlist, [g + h, adj[i]])

            # print('')
            heapify(openlist)
            if (end[0], end[1], end[2]) in info:
                if info[(end[0], end[1], end[2])][0]:
                    break

        path = []
        seg = 0
        via = 0
        if (end[0], end[1], end[2]) in info:
            if info[(end[0], end[1], end[2])][0]:
                path, seg, via = self.update_path(info, start, end,self.demand,1)
        return path, seg, via


if __name__ == '__main__':

    parameter = data_generator.dic_generate(
        32, 32, 2,
        4.0, 4.0,
        [1.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
        10.0,
        10.0,
        {},
        100, 5
    )

    for item in parameter:
        print(item, parameter[item])
    env = Environment(parameter)
    env.load_parameters()
    env.generate_capacity()

    astar = Astar(env)

    total_seg = 0
    total_via = 0

    pinpairlist = env.pinpairs.copy()

    for pinpair in pinpairlist:
        start = pinpair[0][:3]
        end = pinpair[1][:3]
        # print(start, end)
        path, seg, via = astar.search_withcap(start, end)
        total_seg += seg
        total_via += via
    print(total_seg, total_via)

    total_seg = 0
    total_via = 0
    for pinpair in pinpairlist:
        start = pinpair[0][:3]
        end = pinpair[1][:3]
        # print(start,end)
        path, seg, via = astar.search_withcap(start, end)
        total_seg += seg
        total_via += via
        print(seg, via, path)
    print(total_seg, total_via)

    total_seg = 0
    total_via = 0
    astar.reset_cap()
    random.shuffle(pinpairlist)
    for pinpair in pinpairlist:
        start = pinpair[0][:3]
        end = pinpair[1][:3]
        # print(start,end)
        path, seg, via = astar.search_withcap(start, end)
        total_seg += seg
        total_via += via
        print(seg, via, path)
    print(total_seg, total_via)
