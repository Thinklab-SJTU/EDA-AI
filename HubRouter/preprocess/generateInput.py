import numpy as np
from PIL import Image
import os
from grid import Grid
import torch
from torch.nn import functional as F

from utils import read_input_info, get_coordinate
import random
random.seed(333)

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


def gene_input(env, nets, routes_info, grid, bm_name, scale=128, mode='clip', 
                clip_mode='random'):

    if not os.path.exists('../dataset/' + bm_name + '_input_' + str(scale)):
        os.mkdir('../dataset/' + bm_name + '_input_' + str(scale))
    if not os.path.exists('../dataset/' + bm_name + '_' + str(scale)):
        os.mkdir('../dataset/' + bm_name + '_' + str(scale))
        
    line = 0
    n_line = len(routes_info)

    a_len, b_len, c_len, d_len = 0, 0, 0, 0
    print()
    while line < n_line:
        net_id = str(routes_info[line][1])
        net_name = routes_info[line][0]

        condition_img = np.zeros((env.n_x, env.n_y, 3), dtype=np.dtype('uint8'))
        condition_img[:, :, 1] = env.planar_cap[:, :, 0]
        condition_img[:, :, 2] = env.planar_cap[:, :, 1]
        for pin in nets[net_name]['pins']:
            condition_img[pin[0], pin[1], 0] = 255

        route_img = torch.zeros(grid[0], grid[1])
        while routes_info[line+1][0] != '!':
            line += 1
            tile0, tile1 = routes_info[line][0].split('-')
            x1, y1, z1 = tile0[1:-1].split(',')
            x1 = int((int(x1) - grid[2]) / grid[4])
            y1 = int((int(y1) - grid[3]) / grid[5])
            z1 = int(z1)
            x2, y2, z2 = tile1[1:-1].split(',')
            x2 = int((int(x2) - grid[2]) / grid[4])
            y2 = int((int(y2) - grid[3]) / grid[5])
            z2 = int(z2)
            if z1 == z2:
                if x1 == x2:
                    env.planar_cap[x1, y1:y2, 1] -= 1
                    route_img[x1, y1:(y2+1)] = 255
                else:
                    env.planar_cap[x1:x2, y1, 0] -= 1
                    route_img[x1:(x2+1), y1] = 255

        line += 2
        
        min_cap = env.planar_cap[:-1, :-1, :].min()
        print(f'{bm_name}: {line}/{n_line}, a: {a_len}/15000, b: {b_len}/15000, ' + \
               f'c: {c_len}/15000, d: {d_len}/15000, min_cap: {min_cap}', end='\r')
        if min_cap > 0:
            continue

        # only one pixel
        if route_img.sum().item() <= 255:
            continue

        if mode == 'clip':
            if not os.path.exists('../dataset/' + bm_name + '_input_' + str(scale) + '/a/'):
                os.mkdir('../dataset/' + bm_name + '_input_' + str(scale) + '/a/')
            if not os.path.exists('../dataset/' + bm_name + '_' + str(scale) + '/a/'):
                os.mkdir('../dataset/' + bm_name + '_' + str(scale) + '/a/')
            if not os.path.exists('../dataset/' + bm_name + '_input_' + str(scale) + '/b/'):
                os.mkdir('../dataset/' + bm_name + '_input_' + str(scale) + '/b/')
            if not os.path.exists('../dataset/' + bm_name + '_' + str(scale) + '/b/'):
                os.mkdir('../dataset/' + bm_name + '_' + str(scale) + '/b/')
            if not os.path.exists('../dataset/' + bm_name + '_input_' + str(scale) + '/c/'):
                os.mkdir('../dataset/' + bm_name + '_input_' + str(scale) + '/c/')
            if not os.path.exists('../dataset/' + bm_name + '_' + str(scale) + '/c/'):
                os.mkdir('../dataset/' + bm_name + '_' + str(scale) + '/c/')
            if not os.path.exists('../dataset/' + bm_name + '_input_' + str(scale) + '/d/'):
                os.mkdir('../dataset/' + bm_name + '_input_' + str(scale) + '/d/')
            if not os.path.exists('../dataset/' + bm_name + '_' + str(scale) + '/d/'):
                os.mkdir('../dataset/' + bm_name + '_' + str(scale) + '/d/')
                
            route_position = np.where(route_img>1)
            pin_location = np.where(condition_img[:, :, 0]>1)
            pin_num = len(nets[net_name]['pins'])
            pin_wl = (pin_location[0].max() - pin_location[0].min()).item()+1 + (pin_location[1].max() - pin_location[1].min()).item()+1
            if (route_position[0].max() - route_position[0].min() >=scale-2) or \
                (route_position[1].max() - route_position[1].min() >=scale-2):
                # unable to accommodate in (scale,scale) canvas
                continue
            else:
                # decide which fold to be put into
                # a: pin_num<=4 and pin_wl<=16
                # b: pin_num>4 and pin_wl<=16
                # c: pin_num<4 and pin_wl>16
                # d: pin_num>4 and pin_wl>16
                # number of images is up to 15000 in each catagory
                if min(a_len, b_len, c_len, d_len) >= 15000:
                    return 
                
                if pin_num<=4 and pin_wl<=16:
                    letter = 'a'
                    if a_len >= 15000:
                        continue
                    else:
                        a_len += 1
                elif pin_num>4 and pin_wl<=16:
                    letter = 'b'
                    if b_len >= 15000:
                        continue
                    else:
                        b_len += 1
                elif pin_num<=4 and pin_wl>16:
                    letter = 'c'
                    if c_len >= 15000:
                        continue
                    else:
                        c_len += 1
                else:
                    letter = 'd'
                    if d_len >= 15000:
                        continue
                    else:
                        d_len += 1
            
                if clip_mode == 'random':
                    # randomly choose the route position
                    random_x = random.randint(max(route_position[0].max() - scale + 1, 0), min(route_position[0].min(), route_img.size(0)-scale))
                    random_y = random.randint(max(route_position[1].max() - scale + 1, 0), min(route_position[1].min(), route_img.size(1)-scale))
                    route_img = route_img[random_x:random_x+scale, random_y:random_y+scale]
                    condition_img = condition_img[random_x:random_x+scale, random_y:random_y+scale, :]
                elif clip_mode == 'up-left':
                    condition_img = torch.tensor(condition_img)
                    if route_position[0].min() == 0:
                        route_img = torch.nn.functional.pad(route_img, (0, 0, 1, 0), mode='constant', value=0)
                        condition_img = torch.nn.functional.pad(condition_img, (0, 0, 0, 0, 1, 0), mode='constant', value=0)
                    else:
                        route_img = route_img[route_position[0].min()-1: route_position[0].min()-1+scale, :]
                        route_img = F.pad(route_img, (0, 0, 0, max(scale - route_img.shape[0], 0)), mode='constant', value=0)
                        condition_img = condition_img[route_position[0].min()-1: route_position[0].min()-1+scale, :, :]
                        condition_img = torch.nn.functional.pad(condition_img, (0, 0, 0, 0, 0, max(scale - condition_img.shape[0], 0)), mode='constant', value=0)
                    if route_position[1].min() == 0:
                        route_img = torch.nn.functional.pad(route_img, (1, 0, 0, 0), mode='constant', value=0)
                        condition_img = torch.nn.functional.pad(condition_img, (0, 0, 1, 0, 0, 0), mode='constant', value=0)
                    else:
                        route_img = route_img[:, route_position[1].min()-1: route_position[1].min()-1+scale]
                        route_img = torch.nn.functional.pad(route_img, (0, max(scale - route_img.shape[1], 0), 0, 0), mode='constant', value=0)
                        condition_img = condition_img[:, route_position[1].min()-1: route_position[1].min()-1+scale, :]
                        condition_img = torch.nn.functional.pad(condition_img, (0, 0, 0, max(scale - condition_img.shape[1], 0), 0, 0), mode='constant', value=0)

        route_img = route_img.numpy().astype(np.dtype('uint8'))
        if mode == 'clip' and clip_mode == 'up-left':
            condition_img = condition_img.numpy().astype(np.dtype('uint8'))
        
        assert route_img.shape[0] == route_img.shape[1] == condition_img.shape[0] == condition_img.shape[1] == scale
        
        condition_img_file = Image.fromarray(condition_img, mode='RGB')
        condition_img_file.save('../dataset/' + bm_name + '_input_' + str(scale) + '/' + letter + '/' + bm_name + '_' + net_id + '.png')
        route_img_file = Image.fromarray(route_img, mode='L')
        route_img_file.save('../dataset/' + bm_name + '_' + str(scale) + '/' + letter + '/' + bm_name + '_' + str(net_id) + '.png')
    return

if __name__ == '__main__':
    # benchmark = ['adaptec1.capo70.3d.35.50.90.gr', 'adaptec2.mpl60.3d.35.20.100.gr',
                 # 'adaptec3.dragon70.3d.30.50.90.gr', 'adaptec4.aplace60.3d.30.50.90.gr', 
                 # 'adaptec5.mfar50.3d.50.20.100.gr', 'bigblue1.capo60.3d.50.10.100.gr',
                 # 'bigblue2.mpl60.3d.40.60.60.gr', 'bigblue3.aplace70.3d.50.10.90.m8.gr', 
                 # 'bigblue4.fastplace70.3d.80.20.80.gr', 'newblue1.ntup50.3d.30.50.90.gr',
                 # 'newblue2.fastplace90.3d.50.20.100.gr', 'newblue3.kraftwerk80.3d.40.50.90.gr', 
                 # 'newblue4.mpl50.3d.40.10.95.gr', 'newblue5.ntup50.3d.40.10.100.gr', 
                 # 'newblue6.mfar80.3d.60.10.100.gr', 'newblue7.kraftwerk70.3d.80.20.82.m8.gr']
    
    benchmark = ['newblue3.kraftwerk80.3d.40.50.90.gr', 'newblue4.mpl50.3d.40.10.95.gr', 
                 'newblue7.kraftwerk70.3d.80.20.82.m8.gr', 'bigblue4.fastplace70.3d.80.20.80.gr']
    grid = {'adaptec1': [324, 324, 136, 136, 35, 35], 'adaptec2': [424, 424, 216, 216, 35, 35],
            'adaptec3': [774, 779, 21, 43, 30, 30], 'adaptec4': [774, 779, 21, 43, 30, 30],
            'adaptec5': [465, 468, 11, 33, 50, 50], 'bigblue1': [227, 227, 128, 128, 50, 50],
            'bigblue2': [468, 471, 16, 56, 40, 40], 'bigblue3': [555, 557, 11, 51, 50, 50], 
            'bigblue4': [403, 405, 0, 18, 80, 80], 
            'newblue1': [399, 399, 218, 218, 30, 30], 'newblue2': [557, 463, 7, 7, 50, 50],
            'newblue3': [973, 1256, 16, 19, 40, 40], 'newblue4': [455, 458, 4, 12, 40, 40], 
            'newblue5': [637, 640, 16, 47, 40, 40], 'newblue6': [463, 464, 6, 46, 60, 60], 
            'newblue7': [488, 490, 0, 9, 80, 80]}
    
    for i in range(len(benchmark)):
        bm_name = benchmark[i].split('.')[0]
        grid_info = read_input_info('./benchmark/' + benchmark[i])
        
        grid_param, nets = get_input_parameters(grid_info)
        env = Grid(grid_param)
        routes_info = read_input_info('./output/' + bm_name)
        gene_input(env, nets, routes_info, grid[bm_name], bm_name, scale=64, mode='clip', clip_mode='random')

