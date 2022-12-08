import numpy as np


class Grid(object):
    def __init__(self, grid_parameters):
        self.n_x = grid_parameters['grid_size'][0]  # #row
        self.n_y = grid_parameters['grid_size'][1]  # #col
        self.vtc_cap = grid_parameters['vertical_capacity']  # vertical capacity
        self.hzt_cap = grid_parameters['horizontal_capacity']  # horizontal capacity
        # self.twl = 0  # total wire length
        # self.tof = 0  # total overflow
        self.grid = self.init_capacity()

    def init_capacity(self):
        capacity = np.zeros((3, self.n_x, self.n_y))
        capacity[1, :, :] = self.hzt_cap
        capacity[2, :, :] = self.vtc_cap
        capacity[1, -1, :] = 0
        capacity[2, :, -1] = 0
        return capacity

# if __name__ == '__main__':
#     filename = '.\\data\\adaptec1.capo70.3d.35.50.90.gr'
#     grid_info = utils.read_input_info(filename)
#     grid_parameters, nets = init.get_input_parameters(grid_info)
#     grid = Grid(grid_parameters)
#     print(len(nets), grid.planar_cap)
