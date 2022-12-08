from .grid import Grid


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


# Parsing input data

def preprocess(grid_info):
    grid_parameters = {'grid_size': [int(grid_info[0][1]), int(grid_info[0][2])],
                       'vertical_capacity': int(grid_info[1][2]),
                       'horizontal_capacity': int(grid_info[2][2])
                       }
    num_net = int(grid_info[3][2])

    line = 4
    nets = []
    for _ in range(num_net):
        if line >= len(grid_info):
            break
        net = {'net_name': str(grid_info[line][0]), 'net_ID': int(grid_info[line][1]),
               'num_pin': int(grid_info[line][2])}

        # skip the nets with more than 1000 pins or all the pins falling in the same tile
        # if net['num_pin'] > 1000:
        #     line += net['num_pin'] + 1
        #     continue

        line += 1
        pins = []
        for __ in range(net['num_pin']):
            pin = [int(grid_info[line][i]) for i in range(2)]
            # here, the type of layer no. is float
            # pins.append([float(grid_info[line][i]) for i in range(3)])
            pins.append(pin)
            line += 1

        is_in_the_same_tile = True
        x1, y1 = pins[0][0], pins[0][1]
        for i in range(1, len(pins)):
            x2, y2 = pins[i][0], pins[i][1]
            if x1 != x2 or y1 != y2:
                is_in_the_same_tile = False
                break
        if is_in_the_same_tile:
            continue

        net['pins'] = pins
        nets.append(pins)
        # nets.append(net)
        # nets[net['net_name']] = net

    return Grid(grid_parameters), nets

