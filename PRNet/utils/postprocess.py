from concurrent.futures import process
import torch
import numpy as np
import collections
from PIL import Image
import time
from ctypes import CDLL, c_int, c_double, byref, c_float, c_bool

so = CDLL("./utils/new_c_utils.so")


# def dfs(grid, r, c, num, m, n, pin, record):
#     flag = False
#     over = False
#     grid[r][c] = num
#     # m, n = len(grid), len(grid[0])
#     if (r * m + c).item() in pin:
#         flag = True
#     record.discard((r * m + c).item())
#     if len(record) == 0:
#         return True, True

#     for x, y in [(r - 1, c), (r + 1, c), (r, c - 1), (r, c + 1)]:
#         if 0 <= x < m and 0 <= y < n and grid[x][y] == 1:
#             temp_flag, temp_over = dfs(grid, x, y, num, m, n, pin, record)
#             flag = flag or temp_flag
#             over = over or temp_over
#             if over:
#                 return True, True
#     if (r * m + c).item() not in pin and not flag:
#         grid[r][c] = 0
#         return False, over
#     return True, over


def postprocess(tensor_A, tensor_B):
    """
    input:
        A: input matrix 3 x 64 x 64
        B: generated route 1 x 64 x 64
    output:
        flag: Boolean -> succeed or not
        grid: success -> output the postprocessed route 1 x 64 x 64
              failue -> False
    """

    bs = tensor_A.shape[0]
    m, n = tensor_A.shape[2], tensor_A.shape[3]
    res_flag = []
    res_grid = []
    res_skipwl = 0

    time1, time2, time3 = 0, 0, 0

    for i in range(bs):

        # start_time = time.time()

        # grid = tensor_B[i, 0].clone()
        grid = tensor_B[i, 0]

        grid = grid.reshape(m * n).astype(np.float64)
        A = tensor_A[i, 0].reshape(m * n).astype(np.float64)

        # time1 += time.time() - start_time
        # start_time = time.time()

        # grid = (c_double * (n*n))(*grid)
        # A = (c_double * (n*n))(*A)
        grid = (c_double * len(grid)).from_buffer(grid)
        A = (c_double * len(A)).from_buffer(A)
        num = c_int(2)
        over = c_bool(False)

        # time2 += time.time() - start_time
        # start_time = time.time()

        skip_wl = so.process(grid, A, m, n, byref(num), byref(over))

        grid = np.array(grid).reshape(m, n)

        # time3 += time.time() - start_time
        # start_time = time.time()

        print(num, skip_wl)
        if num == 2:
            res_flag.append(True)
            res_grid.append(grid)
            res_skipwl += 0
            continue

        res_flag.append(False)
        res_grid.append(grid)
        res_skipwl += skip_wl

    # print(1, time1)
    # print(2, time2)
    # print(3, time3)
    # print("\n")

    return res_skipwl, res_grid


def proc(tensor_A, tensor_B):
    # start_time = time.time()
    skip_wl, grids = postprocess(tensor_A, tensor_B)
    gs = np.array(grids)
    wire_length = np.sum(gs) - gs.shape[0]
    # print(4, time.time() - start_time)

    # return skips, np.sum([grids[i] for i in range(len(grids))]) - len(grids), grids
    return skip_wl, wire_length, grids
