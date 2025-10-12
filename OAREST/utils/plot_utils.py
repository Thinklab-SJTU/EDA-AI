import numpy as np
import matplotlib.pyplot as plt
from utils.obstacle_utils import *

# edge_color = 'black'
# edge_width = .5
# term_color = 'black'
# aux_color = 'blue'
# term_size = 4
# steiner_color = 'black'

# def plot_obstacle_rest(input, output, obstacle, ax, degree):
#     x_low, y_low, x_high, y_high = [], [], [], []
#     existed_rec = []
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
#     for i in range(len(obstacle)):
#         x_min, y_min, x_max, y_max = obstacle[i]
#         if points_in_rec(input, obstacle[i]) or rec_overlap(existed_rec, obstacle[i]):
#             continue
#         rec = plt.Rectangle(
#                 (x_min, y_min), # (x,y) left-down corner
#                 x_max - x_min, # width
#                 y_max - y_min, # height
#                 color='maroon',
#                 alpha=0.5)
#         ax.add_patch(rec)
#         existed_rec.append(obstacle[i])
#     for i in range(len(x_low)):
#         plt.plot([x_low[i], x_high[i]], [input[i][1], input[i][1]], '-', color=edge_color, linewidth=edge_width)
#         plt.plot([input[i][0], input[i][0]], [y_low[i], y_high[i]], '-', color=edge_color, linewidth=edge_width)
#     plt.plot(list(input[:degree,0]), list(input[:degree,1]), 's', color=term_color, markerfacecolor='black', markersize=term_size, markeredgewidth=edge_width)
#     plt.plot(list(input[degree:,0]), list(input[degree:,1]), 's', color=aux_color, markerfacecolor=aux_color, markersize=term_size, markeredgewidth=edge_width)
#     plt.xlim(-0.05, 1.05)
#     plt.ylim(-0.05, 1.05)
#     # fig.savefig('rsmt.pdf')
#     return ax

# def plot_gst_rsmt(terms, sps, edges):
#     degree = len(terms)
#     if len(sps) != 0:
#         points = np.concatenate([terms, sps], 0)
#     else:
#         points = terms
#     for i in range(len(edges)):
#         u = edges[i][0]
#         v = edges[i][1]
#         plt.plot([points[u][0], points[u][0]], [points[u][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
#         plt.plot([points[u][0], points[v][0]], [points[v][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
#     plt.plot([terms[i][0] for i in range(degree)], [terms[i][1] for i in range(degree)], 's', markerfacecolor='black', color=term_color, markersize=term_size, markeredgewidth=edge_width)
#     plt.xlim(-0.05, 1.05)
#     plt.ylim(-0.05, 1.05)      
    
# def plot_obstacle_gst_rsmt(terms, sps, edges, obstacle, ax):
#     degree = len(terms)
#     if len(sps) != 0:
#         points = np.concatenate([terms, sps], 0)
#     else:
#         points = terms
#     existed_rec = []
#     for i in range(len(edges)):
#         u = edges[i][0]
#         v = edges[i][1]
#         plt.plot([points[u][0], points[u][0]], [points[u][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
#         plt.plot([points[u][0], points[v][0]], [points[v][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
#     for i in range(len(obstacle)):
#         x_min, y_min, x_max, y_max = obstacle[i]
#         if points_in_rec(terms, obstacle[i]) or rec_overlap(existed_rec, obstacle[i]):
#             continue
#         rec = plt.Rectangle(
#                 (x_min, y_min), # (x,y) left-down corner
#                 x_max - x_min, # width
#                 y_max - y_min, # height
#                 color='maroon',
#                 alpha=0.5)
#         ax.add_patch(rec)
#         existed_rec.append(obstacle[i])
#     plt.plot([terms[i][0] for i in range(degree)], [terms[i][1] for i in range(degree)], 's', markerfacecolor='black', color=term_color, markersize=term_size, markeredgewidth=edge_width)
#     plt.xlim(-0.05, 1.05)
#     plt.ylim(-0.05, 1.05)  

import numpy as np
import matplotlib.pyplot as plt
from utils.obstacle_utils import *

# Color and style settings
# edge_color = 'darkblue'
# term_color = 'red'
# aux_color = 'green'
# term_size = 6
# steiner_color = 'orange'
# obstacle_color = 'maroon'
# obstacle_alpha = 0.6

# edge_color = '#4E79A7'
# edge_width = 1.5
# term_color = '#F28E2B'
# aux_color = '#76B7B2'
# term_size = 5
# obstacle_color = '#E15759'
# obstacle_alpha = 0.7
# steiner_color = '#59A14F'

# edge_color = '#4E79A7'
edge_color = '#606060'
edge_width = 1.5
# term_color = '#F28E2B'
term_color = '#FFA500'
# aux_color = '#76B7B2'
aux_color = '#87CEEB'
term_size = 5
# obstacle_color = '#E15759'
obstacle_color = '#A9A9A9'
obstacle_alpha = 0.7
# obstacle_color = '#A9A9A9'
# obstacle_alpha = 0.7

def plot_obstacle_rest(input, output, obstacle, ax, degree):
    x_low, y_low, x_high, y_high = [], [], [], []
    existed_rec = []

    # Calculate low and high bounds
    for i in range(len(input)):
        x_low.append(input[i][0])
        y_low.append(input[i][1])
        x_high.append(input[i][0])
        y_high.append(input[i][1])
    
    for i in range(int(len(output) / 2)):
        x_idx = output[2 * i]
        y_idx = output[2 * i + 1]
        x = input[x_idx][0]
        y = input[y_idx][1]
        y_low[x_idx] = min(y_low[x_idx], y)
        y_high[x_idx] = max(y_high[x_idx], y)
        x_low[y_idx] = min(x_low[y_idx], x)
        x_high[y_idx] = max(x_high[y_idx], x)
    
    # Plot obstacles
    for i in range(len(obstacle)):
        x_min, y_min, x_max, y_max = obstacle[i]
        if points_in_rec(input, obstacle[i]) or rec_overlap(existed_rec, obstacle[i]):
            continue
        rec = plt.Rectangle(
            (x_min, y_min),  # (x,y) left-down corner
            x_max - x_min,  # width
            y_max - y_min,  # height
            color=obstacle_color,
            alpha=obstacle_alpha)
        ax.add_patch(rec)
        existed_rec.append(obstacle[i])
    
    # Plot the edges
    for i in range(len(x_low)):
        plt.plot([x_low[i], x_high[i]], [input[i][1], input[i][1]], '-', color=edge_color, linewidth=edge_width)
        plt.plot([input[i][0], input[i][0]], [y_low[i], y_high[i]], '-', color=edge_color, linewidth=edge_width)
    
    # Plot terminals and auxiliary points
    plt.plot(list(input[:degree, 0]), list(input[:degree, 1]), 's', color=term_color, markerfacecolor=term_color, markersize=term_size, markeredgewidth=edge_width)
    plt.plot(list(input[degree:, 0]), list(input[degree:, 1]), 'o', color=aux_color, markerfacecolor=aux_color, markersize=term_size, markeredgewidth=edge_width)
    
    # Set limits and grid
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    # plt.grid(True, linestyle='--', alpha=0.7)

    return ax

# def plot_gst_rsmt(terms, sps, edges):
#     degree = len(terms)
#     if len(sps) != 0:
#         points = np.concatenate([terms, sps], 0)
#     else:
#         points = terms
    
#     # Plot edges
#     for i in range(len(edges)):
#         u = edges[i][0]
#         v = edges[i][1]
#         plt.plot([points[u][0], points[u][0]], [points[u][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
#         plt.plot([points[u][0], points[v][0]], [points[v][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
    
#     # Plot terminals
#     plt.plot([terms[i][0] for i in range(degree)], [terms[i][1] for i in range(degree)], 's', markerfacecolor=term_color, color=term_color, markersize=term_size, markeredgewidth=edge_width)
    
#     # Set limits and grid
#     plt.xlim(-0.05, 1.05)
#     plt.ylim(-0.05, 1.05)
#     # plt.grid(True, linestyle='--', alpha=0.7)

def plot_obstacle_gst_rsmt(terms, sps, edges, obstacle, ax):
    degree = len(terms)
    if len(sps) != 0:
        points = np.concatenate([terms, sps], 0)
    else:
        points = terms
    existed_rec = []
    
    # Plot edges
    for i in range(len(edges)):
        u = edges[i][0]
        v = edges[i][1]
        plt.plot([points[u][0], points[u][0]], [points[u][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
        plt.plot([points[u][0], points[v][0]], [points[v][1], points[v][1]], '-', color=edge_color, linewidth=edge_width)
    
    # Plot obstacles
    for i in range(len(obstacle)):
        x_min, y_min, x_max, y_max = obstacle[i]
        if points_in_rec(terms, obstacle[i]) or rec_overlap(existed_rec, obstacle[i]):
            continue
        rec = plt.Rectangle(
            (x_min, y_min),  # (x,y) left-down corner
            x_max - x_min,  # width
            y_max - y_min,  # height
            color=obstacle_color,
            alpha=obstacle_alpha)
        ax.add_patch(rec)
        existed_rec.append(obstacle[i])
    
    # Plot terminals
    plt.plot([terms[i][0] for i in range(degree)], [terms[i][1] for i in range(degree)], 's', markerfacecolor=term_color, color=term_color, markersize=term_size, markeredgewidth=edge_width)
    
    # Set limits and grid
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    # plt.grid(True, linestyle='--', alpha=0.7)

    return ax
