import torch 
import numpy as np

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

def get_coordinate(position, origin, tile_width, tile_height):
    x = int((position[0] - origin[0]) / tile_width)
    y = int((position[1] - origin[1]) / tile_height)
    z = int(position[2])
    return x, y, z

def get_coordinate_planar(position, origin, tile_width, tile_height):
    x = int((position[0] - origin[0]) / tile_width)
    y = int((position[1] - origin[1]) / tile_height)
    return x, y

def get_input_parameters(grid_info):
    n_layer = int(grid_info[0][3])
    grid_parameters = {'grid_size': [int(grid_info[0][1]), int(grid_info[0][2]), n_layer],
                       'vertical_capacity': [int(grid_info[1][i + 2]) for i in range(n_layer)],
                       'horizontal_capacity': [int(grid_info[2][i + 2]) for i in range(n_layer)],
                       'min_width': [int(grid_info[3][i + 2]) for i in range(n_layer)],
                       'min_spacing': [int(grid_info[4][i + 2]) for i in range(n_layer)],
                       'via_spacing': [int(grid_info[5][i + 2]) for i in range(n_layer)],
                       'origin': [int(grid_info[6][0]), int(grid_info[6][1])],
                       'tile_width': int(grid_info[6][2]),
                       'tile_height': int(grid_info[6][3])}

    num_net = int(grid_info[7][2])

    line = 8
    nets = {}
    for _ in range(num_net):
        try:
            net = {'net_name': str(grid_info[line][0]), 'net_ID': int(grid_info[line][1]),
               'num_pin': int(grid_info[line][2]), 'min_width': int(grid_info[line][3])}
        except:
            net = {'net_name': str(grid_info[line][0]), 'net_ID': int(grid_info[line][1]),
               'num_pin': int(grid_info[line][2])}

        line += 1
        pins = []
        for __ in range(net['num_pin']):
            pin = [int(grid_info[line][i]) for i in range(2)]
            pin.append(int(grid_info[line][2]))

            pin = get_coordinate(pin, grid_parameters['origin'], grid_parameters['tile_width'],
                                 grid_parameters['tile_height'])
            pins.append(pin)
            line += 1

        is_in_the_same_tile = True
        x1, y1, z1 = pins[0][0], pins[0][1], pins[0][2]
        for i in range(len(pins)):
            x2, y2, z2 = pins[i][0], pins[i][1], pins[i][2]
            if x1 != x2 or y1 != y2 or z1 != z2:
                is_in_the_same_tile = False
                break
        if is_in_the_same_tile:
            continue

        net['pins'] = pins
        nets[net['net_name']] = net

    num_capacity_adjustment = int(grid_info[line][0])
    line += 1

    reduced_capacity = []
    for _ in range(num_capacity_adjustment):
        reduced_capacity.append([int(grid_info[line][i]) for i in range(7)])
        line += 1

    grid_parameters['reduced_capacity'] = reduced_capacity

    return grid_parameters, nets

def get_input_parameters_ibm(grid_info):
    n_layer = int(grid_info[0][3])
    grid_parameters = {'grid_size': [int(grid_info[0][1]), int(grid_info[0][2]), n_layer],
                       'vertical_capacity': [int(grid_info[1][i + 2]) for i in range(n_layer)],
                       'horizontal_capacity': [int(grid_info[2][i + 2]) for i in range(n_layer)],
                       'min_width': [int(grid_info[3][i + 2]) for i in range(n_layer)],
                       'min_spacing': [int(grid_info[4][i + 2]) for i in range(n_layer)],
                       'via_spacing': [int(grid_info[5][i + 2]) for i in range(n_layer)],
                       'origin': [int(grid_info[6][0]), int(grid_info[6][1])],
                       'tile_width': int(grid_info[6][2]),
                       'tile_height': int(grid_info[6][3])}

    num_net = int(grid_info[7][2])

    line = 8
    nets = {}
    for _ in range(num_net):
        net = {'net_name': str(grid_info[line][0]), 'net_ID': int(grid_info[line][1]),
               'num_pin': int(grid_info[line][2])}

        line += 1
        pins = []
        for _ in range(net['num_pin']):
            pin = [int(grid_info[line][i]) for i in range(2)]

            pin = get_coordinate_planar(pin, grid_parameters['origin'], grid_parameters['tile_width'],
                                        grid_parameters['tile_height'])
            pins.append(pin)
            line += 1

        is_in_the_same_tile = True
        x1, y1 = pins[0][0], pins[0][1]
        for i in range(len(pins)):
            x2, y2 = pins[i][0], pins[i][1]
            if x1 != x2 or y1 != y2:
                is_in_the_same_tile = False
                break
        # if is_in_the_same_tile:
            # continue

        net['pins'] = pins
        net['same_tile'] = is_in_the_same_tile
        nets[net['net_name']] = net
        
    try:
        num_capacity_adjustment = int(grid_info[line][0])
        line += 1

        reduced_capacity = []
        for _ in range(num_capacity_adjustment):
            reduced_capacity.append([int(grid_info[line][i]) for i in range(7)])
            line += 1

        grid_parameters['reduced_capacity'] = reduced_capacity
    except:
        reduced_capacity = []
        grid_parameters['reduced_capacity'] = reduced_capacity

    return grid_parameters, nets