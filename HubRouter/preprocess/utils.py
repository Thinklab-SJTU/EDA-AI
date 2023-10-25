import numpy as np
import os


def read_input_info(filename):
    file = open(filename, 'r')
    grid_info = {}
    i = 0
    for line in file:
        if not line.strip():
            continue
        else:
            grid_info[i] = line.split()
        i += 1
    file.close()
    return grid_info


def is_different(pos1, pos2):
    if len(pos1) != len(pos2):
        return True
    for i in range(len(pos1)):
        if pos1[i] != pos2[i]:
            return True
    return False


def compare(tp1, tp2):
    measure1 = np.abs(tp1['pin_pair'][0][0] - tp1['pin_pair'][1][0]) \
               + np.abs(tp1['pin_pair'][0][1] - tp1['pin_pair'][1][1]) \
               + np.abs(tp1['pin_pair'][0][2] - tp1['pin_pair'][1][2])
    measure2 = np.abs(tp2['pin_pair'][0][0] - tp2['pin_pair'][1][0]) \
               + np.abs(tp2['pin_pair'][0][1] - tp2['pin_pair'][1][1]) \
               + + np.abs(tp2['pin_pair'][0][2] - tp2['pin_pair'][1][2])
    return measure1 - measure2


def compare_tile(tp1, tp2):
    measure1 = np.abs(tp1['pin_pair_tile'][0][0] - tp1['pin_pair_tile'][1][0]) \
               + np.abs(tp1['pin_pair_tile'][0][1] - tp1['pin_pair_tile'][1][1]) \
               + np.abs(tp1['pin_pair_tile'][0][2] - tp1['pin_pair_tile'][1][2])
    measure2 = np.abs(tp2['pin_pair_tile'][0][0] - tp2['pin_pair_tile'][1][0]) \
               + np.abs(tp2['pin_pair_tile'][0][1] - tp2['pin_pair_tile'][1][1]) \
               + + np.abs(tp2['pin_pair_tile'][0][2] - tp2['pin_pair_tile'][1][2])
    return measure1 - measure2


def slice_net_min_width(nets):
    min_width_dict = {}
    for net in nets:
        min_width_dict[nets[net]['net_ID']] = nets[net]['min_width']
    return min_width_dict


def cal_ovfl(grid):
    tof = 0.
    max_of = 0.
    for i in range(grid.n_x):
        for j in range(grid.n_y):
            for k in range(grid.n_layer):
                for m in range(6):
                    if grid.capacity[i, j, k, m] < 0:
                        tof += -grid.capacity[i, j, k, m]
                        if -grid.capacity[i, j, j, m] > max_of:
                            max_of = -grid.capacity[i, j, k, m]
    return tof / 2., max_of


def output(result, filename):
    # f = open(filename + '.r', 'a')
    # for net_name in result:
    #     for route in result[net_name]:
    #
    return


def get_coordinate(position, origin, tile_width, tile_height):
    x = int((position[0] - origin[0]) / tile_width)
    y = int((position[1] - origin[1]) / tile_height)
    z = int(position[2])
    return x, y, z

# if __name__ == '__main__':
#
