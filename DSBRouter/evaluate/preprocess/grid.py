import numpy as np
import utils


class Grid(object):
    def __init__(self, grid_parameters):
        self.n_x = grid_parameters['grid_size'][0]  # #row
        self.n_y = grid_parameters['grid_size'][1]  # #col
        self.n_layer = grid_parameters['grid_size'][2]  # #layer
        self.vtc_cap = grid_parameters['vertical_capacity']  # vertical capacity
        self.hzt_cap = grid_parameters['horizontal_capacity']  # horizontal capacity
        self.via_spc = grid_parameters['via_spacing']  # via spacing
        self.min_width = grid_parameters['min_width']  # minimum width of a wire
        self.min_spc = grid_parameters['min_spacing']  # minimum spacing of wires
        self.tile_width = grid_parameters['tile_width']  # tile width
        self.tile_height = grid_parameters['tile_height']  # tile height
        self.origin = grid_parameters['origin']  # origin of each layer
        # self.init_tile = []  # start tile of the two-pin problem
        # self.tmn_tile = []  # end tile of the two-pin problem
        # self.position = []  # current position of the agent
        self.net_width = [self.min_width[i] + self.min_spc[i] for i in range(self.n_layer)]
        self.n_step = 0  # #step the agent has taken
        self.iteration = 0  # iteration this two-pin problem has taken
        self.max_iteration = 0
        self.twl = 0  # total wire length
        self.tof = 0  # total overflow
        self.max_of = 0
        self.via = 0
        self.original_cap = self.init_capacity(grid_parameters['reduced_capacity'])
        # self.cubic_cap = self.original_cap.copy()
        self.original_planar_cap = self.proj_cap()
        self.planar_cap = self.original_planar_cap.copy()

        # print('Grid established.')
        # print('-----------------------------------------------------------')
        return

    def step(self, action, min_width):
        # self.iteration += 1
        # self.n_step += 1
        # x = self.position[0]
        # y = self.position[1]
        # z = self.position[2]
        next_state = self.get_new_state(action, min_width)
        # a = (self.original_cap[x, y, z, action] - self.capacity[x, y, z, action]) / self.net_width[z] + 1.
        # b = self.original_cap[x, y, z, action] / self.net_width[z] + 0.1
        penalty = 100  # (a / b) ** 3
        reward = -1 - penalty
        done = False
        if not utils.is_different(self.position, self.tmn_tile):
            reward += 1000.  # 100., 500., 1000.
            done = True
        elif self.iteration == self.max_iteration:  # 1000, 2000, 10000, increasing with #step(10 * int(np.log10(
            # self.n_step)))
            done = True
        return next_state, reward, done

    def reset(self):
        # self.init_tile = []
        # self.tmn_tile = []
        # self.position = []
        self.n_step = 0
        self.iteration = 0
        self.max_iteration = 0
        self.twl = 0
        self.tof = 0
        self.via = 0
        self.planar_cap = self.original_planar_cap.copy()

    # def set_routing_param(self, pin_pair_tile):
    #     # self.init_tile = pin_pair_tile[0]
    #     # self.tmn_tile = pin_pair_tile[1]
    #     self.max_iteration = 10 * (np.abs(self.init_tile[0] - self.tmn_tile[0])
    #                                + np.abs(self.init_tile[1] - self.tmn_tile[1])
    #                                + np.abs(self.init_tile[2] - self.tmn_tile[2]))

    # def reset_pin_pair(self):
    #     # self.position = []
    #     self.n_step = 0
    #     self.iteration = 0
    #     self.twl = 0
    #     self.tof = 0
    #     self.via = 0
    #     self.capacity = self.original_cap.copy()

    def init_capacity(self, reduced_capacity):
        capacity = np.zeros((self.n_x, self.n_y, self.n_layer, 6))

        # Calculate Available NumNet in each direction
        # vertical_cap_path = []
        # horizontal_cap_path = []
        # for i in range(self.n_layer):
        #     vertical_cap_path.append(self.vtc_cap[i] / (self.min_width[i] + self.min_spc[i]))
        #     horizontal_cap_path.append(self.hzt_cap[i] / (self.min_width[i] + self.min_spc[i]))

        # Apply available NumNet to grid capacity variables
        # for i in range(self.n_layer):
        #     capacity[:, :, i, 0] = capacity[:, :, i, 1] = horizontal_cap_path[i]
        #     capacity[:, :, i, 2] = capacity[:, :, i, 3] = vertical_cap_path[i]

        for i in range(self.n_layer):
            capacity[:, :, i, 0] = capacity[:, :, i, 1] = self.hzt_cap[i] / self.net_width[i]
            capacity[:, :, i, 2] = capacity[:, :, i, 3] = self.vtc_cap[i] / self.net_width[i]

        # Suppose Via Capacity to be large enough
        capacity[:, :, :, 4] = capacity[:, :, :, 5] = 20  # 10, 20, 50, ...

        # Impose capacity adjustments on the initial capacity
        # (action:direction) : (0:+x), (1:-x), (2:+y), (3:-y), (4:+z), (5:-z)
        dirc = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        for r in reduced_capacity:
            tmp_dirc1 = [r[i + 3] - r[i] for i in range(3)]
            tmp_dirc2 = [r[i] - r[i + 3] for i in range(3)]
            idx1 = dirc.index(tmp_dirc1)
            idx2 = dirc.index(tmp_dirc2)
            capacity[r[0], r[1], r[2], idx1] = capacity[r[3], r[4], r[5], idx2] = int(r[6]) / self.net_width[0]

        # Remove capacity outside the grid
        capacity[0, :, :, 1] = capacity[self.n_x - 1, :, :, 0] = 0
        capacity[:, 0, :, 3] = capacity[:, self.n_y - 1, :, 2] = 0
        capacity[:, :, 0, 5] = capacity[:, :, self.n_layer - 1, 4] = 0

        return capacity

    def reset_state(self):
        # self.position = self.init_tile.copy()
        self.iteration = 0
        state_list = []
        # state_list.extend(self.position)  # current position
        state_list.extend([(self.tmn_tile[i] - self.position[i]) for i in range(3)])  # displacement to the end
        state_list.extend(self.capacity[self.position[0], self.position[1], self.position[2], :])
        # state_list.append(self.n_step)
        state = np.array(state_list, dtype=float)
        return state

    def get_new_state(self, action, min_width):
        dirc = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        ex_position = self.position.copy()
        self.position[0] = min(self.n_x - 1, max(0, self.position[0] + dirc[action][0]))
        self.position[1] = min(self.n_y - 1, max(0, self.position[1] + dirc[action][1]))
        self.position[2] = min(self.n_layer - 1, max(0, self.position[2] + dirc[action][2]))

        new_state_list = []
        new_state_list.extend([(self.tmn_tile[i] - self.position[i]) for i in range(3)])  # displacement to the end
        self.update_capacity(action, self.position[0], self.position[1], self.position[2],
                             utils.is_different(ex_position, self.position), min_width)
        new_state_list.extend(self.capacity[self.position[0], self.position[1], self.position[2], :])
        new_state = np.array(new_state_list, dtype=float)
        return new_state

    def update_capacity(self, action, x, y, z, is_different, min_width):
        if is_different:
            if action == 0:
                if self.capacity[x, y, z, 1] <= 0:
                    self.tof += self.min_spc[z] + min_width
                self.capacity[x, y, z, 1] -= self.min_spc[z] + min_width
                self.capacity[x - 1, y, z, 0] -= self.min_spc[z] + min_width
                self.twl += 1
            elif action == 1:
                if self.capacity[x, y, z, 0] <= 0:
                    self.tof += self.min_spc[z] + min_width
                self.capacity[x, y, z, 0] -= self.min_spc[z] + min_width
                self.capacity[x + 1, y, z, 1] -= self.min_spc[z] + min_width
                self.twl += 1
            elif action == 2:
                if self.capacity[x, y, z, 3] <= 0:
                    self.tof += self.min_spc[z] + min_width
                self.capacity[x, y, z, 3] -= self.min_spc[z] + min_width
                self.capacity[x, y - 1, z, 2] -= self.min_spc[z] + min_width
                self.twl += 1
            elif action == 3:
                if self.capacity[x, y, z, 2] <= 0:
                    self.tof += self.min_spc[z] + min_width
                self.capacity[x, y, z, 2] -= self.min_spc[z] + min_width
                self.capacity[x, y + 1, z, 3] -= self.min_spc[z] + min_width
                self.twl += 1
            elif action == 4:
                if self.capacity[x, y, z, 5] <= 0:
                    self.tof += 1.
                self.capacity[x, y, z, 5] -= 1.
                self.capacity[x, y, z - 1, 4] -= 1.
                self.via += 1
                self.twl += 1
            elif action == 5:
                if self.capacity[x, y, z, 4] <= 0:
                    self.tof += 1.
                self.capacity[x, y, z, 4] -= 1.
                self.capacity[x, y, z + 1, 5] -= 1.
                self.via += 1
                self.twl += 1

    def get_coordinate(self, position):
        x = int((position[0] - self.origin[0]) / self.tile_width)
        y = int((position[1] - self.origin[1]) / self.tile_height)
        z = int(position[2])
        return [x, y, z]

    def proj_cap(self):
        planar_capacity = np.zeros((self.n_x, self.n_y, 2))
        for i in range(self.n_x):
            for j in range(self.n_y):
                cap_h = 0
                cap_v = 0
                for k in range(self.n_layer):
                    cap_h += self.original_cap[i][j][k][0]
                    cap_v += self.original_cap[i][j][k][2]
                planar_capacity[i][j][0] = cap_h
                planar_capacity[i][j][1] = cap_v
        return planar_capacity


# if __name__ == '__main__':
#     filename = '.\\data\\adaptec1.capo70.3d.35.50.90.gr'
#     grid_info = utils.read_input_info(filename)
#     grid_parameters, nets = init.get_input_parameters(grid_info)
#     grid = Grid(grid_parameters)
#     print(len(nets), grid.planar_cap)
