import numpy as np
import ctypes
import os
from collections import Counter
from utils import count_overlap

def process_graphs_vectorized(batch_edges, degree):
    batch_size, group_size = batch_edges.shape
    edges_per_graph = group_size // 2
    edges = batch_edges.reshape(batch_size, edges_per_graph, 2)
    edges_mask = edges[...,0] != edges[...,1]  # shape: (batch_size, edges_per_graph)
    results = []
    for i in range(batch_size):
        current_edges = edges[i][edges_mask[i]] 
        while True:
            unique, counts = np.unique(current_edges.flatten(), return_counts=True)
            to_remove = unique[np.logical_and(counts == 1, unique >= degree)]
            if len(to_remove) == 0:
                break
            mask = ~np.any(np.isin(current_edges, to_remove), axis=1)
            new_edges = current_edges[mask]
            
            if new_edges.shape == current_edges.shape:
                break
            current_edges = new_edges
        results.append(current_edges.flatten())
    return results

class ArrayType:
    def __init__(self, type):
        self.type = type

    def from_param(self, param):
        typename = type(param).__name__
        if hasattr(self, 'from_' + typename):
            return getattr(self, 'from_' + typename)(param)
        elif isinstance(param, ctypes.Array):
            return param
        else:
            raise TypeError('Can\'t convert % s' % typename)

    # Cast from lists / tuples
    def from_list(self, param):
        val = ((self.type)*len(param))(*param)
        return val

    from_tuple = from_list
    from_array = from_list
    from_ndarray = from_list

class Evaluator:
    def __init__(self, path):
        eval_mod = ctypes.cdll.LoadLibrary(path)
        DoubleArray = ArrayType(ctypes.c_double)
        IntArray = ArrayType(ctypes.c_int)
        
        eval_mod.gst_open_geosteiner()
        
        eval_func = eval_mod.eval
        eval_func.argtypes = (DoubleArray, IntArray, ctypes.c_int)
        eval_func.restype = ctypes.c_double
        self.eval_func = eval_func
        
        gst_rsmt_func = eval_mod.call_gst_rsmt
        gst_rsmt_func.argtypes = (
            ctypes.c_int, DoubleArray, ctypes.POINTER(ctypes.c_double), 
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double), 
            ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int))
        gst_rsmt_func.restype = ctypes.c_int
        self.gst_rsmt_func = gst_rsmt_func

    def eval_batch(self, input_batch, output_batch, degree):
        lengths = []
        batch_size = len(input_batch)
        for i in range(batch_size):
            lengths.append(self.eval_func(input_batch[i].reshape(-1), output_batch[i], degree))
        return np.array(lengths)


    def eval_obstacle_batch(self, input_batch, input_obstacle_batch, output_batch, degree):
        lengths_oa = []
        overlaps = []
        all_new_input = []
        all_new_output = []

        batch_size = len(input_batch)
        output_batch = process_graphs_vectorized(output_batch, degree)
        for i in range(batch_size):
            test_case = input_batch[i]
            start_index = np.where(test_case == -1)[0]
            start_index = start_index[0] if start_index.size > 0 else degree

            ob_start_index = np.where(input_obstacle_batch[i] == -1)[0]
            ob_start_index = ob_start_index[0] if ob_start_index.size > 0 else degree

            # find the max index of output_list that covers all points
            output_list = output_batch[i].tolist()
            # print(i, ' before: ', output_batch[i])
            # output = output_batch[i].reshape(-1, 2)
            # output = output[output[:,0] != output[:,1]]
            # counts = Counter(output.flatten())
            # to_remove = {k for k, v in counts.items() if v == 1 and k >= start_index}
            # while len(to_remove) != 0:
            #     mask = np.any(np.isin(output, list(to_remove)), axis=1)
            #     output = output[~mask]
            #     counts = Counter(output.flatten())
            #     to_remove = {k for k, v in counts.items() if v == 1 and k >= start_index}
            # output_list = output.flatten().tolist()
            # print(i, ' after: ', output_list)
            covered_indices = [
                output_list.index(k) for k in range(min(degree, start_index))
            ]
            max_index_input = max(covered_indices) if covered_indices else 0
            num_RES = (max_index_input + 2) // 2

            # construct the mapping of obstacles
            tmp_pointindex_to_newindex_dict = {}
            obstacle_points = []
            index = start_index

            for j in range(num_RES * 2):
                out_val = output_list[j]
                if out_val >= degree:
                    if out_val not in tmp_pointindex_to_newindex_dict:
                        obstacle_point_id = out_val - degree
                        obstacle_id = obstacle_point_id // 4
                        obstacle_order = obstacle_point_id % 4
                        obstacle = input_obstacle_batch[i][obstacle_id]

                        # add obstacle points
                        if obstacle_order == 0: # left-down
                            obstacle_points.append([obstacle[0], obstacle[1]])
                        elif obstacle_order == 1: # left-up
                            obstacle_points.append([obstacle[0], obstacle[3]])
                        elif obstacle_order == 2: # right-up
                            obstacle_points.append([obstacle[2], obstacle[3]])
                        elif obstacle_order == 3: # right-down
                            obstacle_points.append([obstacle[2], obstacle[1]])

                        tmp_pointindex_to_newindex_dict[out_val] = index
                        output_list[j] = index
                        index += 1
                    else:
                        output_list[j] = tmp_pointindex_to_newindex_dict[out_val]

            output_list = output_list[:num_RES * 2]

            # construct new inputs
            if obstacle_points:
                obstacle_points = np.array(obstacle_points)
                new_input = np.concatenate(
                    [test_case[:start_index], obstacle_points]
                )
            else:
                new_input = test_case[:start_index]

            all_new_input.append(new_input)
            all_new_output.append(output_list)

            # compute overlap
            overlap = count_overlap(
                new_input, input_obstacle_batch[i][:ob_start_index], output_list
            )

            # compute length
            lengths_oa.append(
                self.eval_func(new_input.reshape(-1), np.array(output_list), index)
            )
            overlaps.append(overlap)

        return np.array(lengths_oa), np.array(overlaps), all_new_input, all_new_output
    
    # def eval_obstacle_batch(self, input_batch, input_obstacle_batch, output_batch, degree):
    #     lengths_oa = []
    #     overlaps = []
    #     all_new_input = []
    #     all_new_output = []

    #     batch_size = len(input_batch)
    #     for i in range(batch_size):
    #         start_index = np.where(input_batch[i] == -1)[0]
    #         start_index = start_index[0] if start_index.size > 0 else degree

    #         ob_start_index = np.where(input_obstacle_batch[i] == -1)[0]
    #         ob_start_index = ob_start_index[0] if ob_start_index.size > 0 else degree

    #         # find the max index of output_list that covers all points
    #         output_list = output_batch[i].tolist()
    #         covered_indices = [
    #             output_list.index(k) for k in range(min(degree, start_index))
    #         ]
    #         max_index_input = max(covered_indices) if covered_indices else 0
    #         num_RES = (max_index_input + 2) // 2

    #         # construct the mapping of obstacles
    #         tmp_pointindex_to_newindex_dict = {}
    #         obstacle_points = []
    #         index = start_index

    #         for j in range(num_RES * 2):
    #             out_val = output_list[j]
    #             if out_val >= degree:
    #                 if out_val not in tmp_pointindex_to_newindex_dict:
    #                     obstacle_point_id = out_val - degree
    #                     obstacle_id = obstacle_point_id // 2
    #                     obstacle_order = obstacle_point_id % 2
    #                     obstacle = input_obstacle_batch[i][obstacle_id]

    #                     # add obstacle points
    #                     if obstacle_order == 0:
    #                         obstacle_points.append([obstacle[0], obstacle[1]])
    #                     else:
    #                         obstacle_points.append([obstacle[2], obstacle[3]])

    #                     tmp_pointindex_to_newindex_dict[out_val] = index
    #                     output_list[j] = index
    #                     index += 1
    #                 else:
    #                     output_list[j] = tmp_pointindex_to_newindex_dict[out_val]

    #         output_list = output_list[:num_RES * 2]

    #         # construct new inputs
    #         if obstacle_points:
    #             obstacle_points = np.array(obstacle_points)
    #             new_input = np.concatenate(
    #                 [input_batch[i][:start_index], obstacle_points]
    #             )
    #         else:
    #             new_input = input_batch[i][:start_index]

    #         all_new_input.append(new_input)
    #         all_new_output.append(output_list)

    #         # compute overlap
    #         overlap = count_overlap(
    #             new_input, input_obstacle_batch[i][:ob_start_index], output_list
    #         )

    #         # compute length
    #         lengths_oa.append(
    #             self.eval_func(new_input.reshape(-1), np.array(output_list), index)
    #         )
    #         overlaps.append(overlap)

    #     return np.array(lengths_oa), np.array(overlaps), all_new_input, all_new_output

    def gst_rsmt(self, inputs):
        degree = len(inputs)
        terms = inputs.flatten()
        length = ctypes.c_double()
        nsps = ctypes.c_int()
        sps = (ctypes.c_double * (degree * 2))()
        sps = ctypes.cast(sps, ctypes.POINTER(ctypes.c_double))
        nedges = ctypes.c_int()
        edges = (ctypes.c_int * (degree * 4))()
        ctypes.cast(edges, ctypes.POINTER(ctypes.c_int))
        self.gst_rsmt_func(degree, terms, ctypes.byref(length), 
            ctypes.byref(nsps), sps, ctypes.byref(nedges), edges)
        sp_list = []
        for i in range(nsps.value):
            sp_list.append([sps[2 * i], sps[2 * i + 1]])
        edge_list = []
        for i in range(nedges.value):
            edge_list.append([edges[2 * i], edges[2 * i + 1]])
        return length.value, sp_list, edge_list

    def gst_obstacle_rsmt(self, inputs, obstacles):
        degree = len(inputs)
        terms = inputs.flatten()
        length = ctypes.c_double()
        nsps = ctypes.c_int()
        sps = (ctypes.c_double * (degree * 2))()
        sps = ctypes.cast(sps, ctypes.POINTER(ctypes.c_double))
        nedges = ctypes.c_int()
        edges = (ctypes.c_int * (degree * 4))()
        ctypes.cast(edges, ctypes.POINTER(ctypes.c_int))
        self.gst_rsmt_func(degree, terms, ctypes.byref(length), 
            ctypes.byref(nsps), sps, ctypes.byref(nedges), edges)
        sp_list = []
        for i in range(nsps.value):
            sp_list.append([sps[2 * i], sps[2 * i + 1]])
        edge_list = []
        for i in range(nedges.value):
            edge_list.append([edges[2 * i], edges[2 * i + 1]])

        # overlaps
        if len(sp_list) != 0:
            inputs = np.concatenate([inputs, np.array(sp_list)], axis=0)
        else:
            pass
        overlaps = count_overlap(inputs, obstacles, np.array(edge_list).flatten())
        return length.value, overlaps, sp_list, edge_list
    
# def horizon_line_in_obstacles(line, y, input_obstacles):
#     overlap = 0
#     for i in range(len(input_obstacles)):
#         if input_obstacles[i][1] < y < input_obstacles[i][3] and \
#             (not (line[0] >= input_obstacles[i][2] or line[1] <= input_obstacles[i][0])):
#             overlap += 1
#     return overlap

# def verticle_line_in_obstacles(line, x, input_obstacles):
#     overlap = 0
#     for i in range(len(input_obstacles)):
#         if input_obstacles[i][0] < x < input_obstacles[i][2] and \
#             (not (line[0] >= input_obstacles[i][3] or line[1] <= input_obstacles[i][1])):
#             overlap += 1
#     return overlap

# def count_overlap(input, input_obstacles, output):
#     # input = new_input
#     # input_obstacles = input_obstacle_batch[i]
#     # output = output_list
#     x_low, y_low, x_high, y_high = [], [], [], []
#     overlap = 0
#     for i in range(len(input)):
#         x_low.append(input[i][0])
#         y_low.append(input[i][1])
#         x_high.append(input[i][0])
#         y_high.append(input[i][1])
#     for i in range(int(len(output) / 2)):
#         x_idx = output[2*i]
#         y_idx = output[2*i+1]
#         x = input[x_idx][0]
#         y = input[y_idx][1]
#         y_low[x_idx] = min(y_low[x_idx], y)
#         y_high[x_idx] = max(y_high[x_idx], y)
#         x_low[y_idx] = min(x_low[y_idx], x)
#         x_high[y_idx] = max(x_high[y_idx], x)
#     for i in range(len(x_low)):
#         if x_low[i] != x_high[i]:
#             line = [x_low[i], x_high[i]]
#             overlap += horizon_line_in_obstacles(line, input[i][1], input_obstacles)
#         if y_low[i] != y_high[i]:
#             line = [y_low[i], y_high[i]]
#             overlap += verticle_line_in_obstacles(line, input[i][0], input_obstacles)
#     return overlap

def plot_rest(input, output):
    x_low, y_low, x_high, y_high = [], [], [], []
    for i in range(len(input)):
        x_low.append(input[i][0])
        y_low.append(input[i][1])
        x_high.append(input[i][0])
        y_high.append(input[i][1])
    for i in range(int(len(output) / 2)):
        x_idx = output[2*i]
        y_idx = output[2*i+1]
        x = input[x_idx][0]
        y = input[y_idx][1]
        y_low[x_idx] = min(y_low[x_idx], y)
        y_high[x_idx] = max(y_high[x_idx], y)
        x_low[y_idx] = min(x_low[y_idx], x)
        x_high[y_idx] = max(x_high[y_idx], x)
    for i in range(len(x_low)):
        plt.plot([x_low[i], x_high[i]], [input[i][1], input[i][1]], '-', color=edge_color, linewidth=edge_width)
        plt.plot([input[i][0], input[i][0]], [y_low[i], y_high[i]], '-', color=edge_color, linewidth=edge_width)
    plt.plot(list(input[:,0]), list(input[:,1]), 's', color=term_color, markerfacecolor='black', markersize=term_size, markeredgewidth=edge_width)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    # fig.savefig('rsmt.pdf')

def plot_obstacle_rest(input, output, obstacle, ax, degree):
    x_low, y_low, x_high, y_high = [], [], [], []
    existed_rec = []
    for i in range(len(input)):
        x_low.append(input[i][0])
        y_low.append(input[i][1])
        x_high.append(input[i][0])
        y_high.append(input[i][1])
    for i in range(int(len(output) / 2)):
        x_idx = output[2*i]
        y_idx = output[2*i+1]
        x = input[x_idx][0]
        y = input[y_idx][1]
        y_low[x_idx] = min(y_low[x_idx], y)
        y_high[x_idx] = max(y_high[x_idx], y)
        x_low[y_idx] = min(x_low[y_idx], x)
        x_high[y_idx] = max(x_high[y_idx], x)
    for i in range(len(obstacle)):
        x_min, y_min, x_max, y_max = obstacle[i]
        if points_in_rec(input, obstacle[i]) or rec_overlap(existed_rec, obstacle[i]):
            continue
        rec = plt.Rectangle(
                (x_min, y_min), # (x,y)矩形左下角
                x_max - x_min, # width长
                y_max - y_min, # height宽
                color='maroon',
                alpha=0.5)
        ax.add_patch(rec)
        existed_rec.append(obstacle[i])
    for i in range(len(x_low)):
        plt.plot([x_low[i], x_high[i]], [input[i][1], input[i][1]], '-', color=edge_color, linewidth=edge_width)
        plt.plot([input[i][0], input[i][0]], [y_low[i], y_high[i]], '-', color=edge_color, linewidth=edge_width)
    plt.plot(list(input[:degree,0]), list(input[:degree,1]), 's', color=term_color, markerfacecolor='black', markersize=term_size, markeredgewidth=edge_width)
    plt.plot(list(input[degree:,0]), list(input[degree:,1]), 's', color=aux_color, markerfacecolor=aux_color, markersize=term_size, markeredgewidth=edge_width)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    # fig.savefig('rsmt.pdf')
    return ax

def plot_gst_rsmt(terms, sps, edges):
    degree = len(terms)
    if len(sps) != 0:
        points = np.concatenate([terms, sps], 0)
    else:
        points = terms
    for i in range(len(edges)):
        u = edges[i][0]
        v = edges[i][1]
        plt.plot([points[u][0], points[u][0]], [points[u][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
        plt.plot([points[u][0], points[v][0]], [points[v][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
    plt.plot([terms[i][0] for i in range(degree)], [terms[i][1] for i in range(degree)], 's', markerfacecolor='black', color=term_color, markersize=term_size, markeredgewidth=edge_width)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)      
    
def plot_obstacle_gst_rsmt(terms, sps, edges, obstacle, ax):
    degree = len(terms)
    if len(sps) != 0:
        points = np.concatenate([terms, sps], 0)
    else:
        points = terms
    existed_rec = []
    for i in range(len(edges)):
        u = edges[i][0]
        v = edges[i][1]
        plt.plot([points[u][0], points[u][0]], [points[u][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
        plt.plot([points[u][0], points[v][0]], [points[v][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
    for i in range(len(obstacle)):
        x_min, y_min, x_max, y_max = obstacle[i]
        if points_in_rec(terms, obstacle[i]) or rec_overlap(existed_rec, obstacle[i]):
            continue
        rec = plt.Rectangle(
                (x_min, y_min), # (x,y)矩形左下角
                x_max - x_min, # width长
                y_max - y_min, # height宽
                color='maroon',
                alpha=0.5)
        ax.add_patch(rec)
        existed_rec.append(obstacle[i])
    plt.plot([terms[i][0] for i in range(degree)], [terms[i][1] for i in range(degree)], 's', markerfacecolor='black', color=term_color, markersize=term_size, markeredgewidth=edge_width)
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)       
    
def save_data(test_data, test_file):
    with open(test_file, 'w') as f:
        for data in test_data:
            f.write(' '.join(['{:.8f} {:.8f}'.format(term[0], term[1]) for term in data]))
            f.write('\n')
            
def read_data(test_file):
    with open(test_file, 'r') as f:
        test_data = []
        lines = f.readlines()
        for line in lines:
            line = line.strip().split(' ')
            data = [float(coord) for coord in line]
            data = np.array(data).reshape([-1, 2])
            test_data.append(data)
    return np.array(test_data)
        

# def transform_inputs(inputs, t):
#     # 0 <= t <= 7
#     xs = inputs[:,:,0]
#     ys = inputs[:,:,1]
#     if t >= 4:
#         temp = xs
#         xs = ys
#         ys = temp
#     if t % 2 == 1:
#         xs = 1 - xs
#     if t % 4 >= 2:
#         ys = 1 - ys
#     xs[np.isclose(xs, 2.0)] = -1
#     ys[np.isclose(ys, 2.0)] = -1
#     return np.stack([xs, ys], -1)
    
# def transform_obstacles(obstacles, t):
#     # 0 <= t <= 7
#     xs1 = obstacles[:,:,0]
#     ys1 = obstacles[:,:,1]
#     xs2 = obstacles[:,:,2]
#     ys2 = obstacles[:,:,3]
#     if t >= 4:
#         xs1, ys1 = ys1, xs1
#         xs2, ys2 = ys2, xs2
#     if t % 2 == 1:
#         xs1, xs2 = 1 - xs2, 1 - xs1
#     if t % 4 >= 2:
#         ys1, ys2 = 1 - ys2, 1 - ys1
#     xs1[np.isclose(xs1, 2.0)] = -1
#     xs2[np.isclose(xs2, 2.0)] = -1
#     ys1[np.isclose(ys1, 2.0)] = -1
#     ys2[np.isclose(ys2, 2.0)] = -1
#     return np.stack([xs1, ys1, xs2, ys2], -1)