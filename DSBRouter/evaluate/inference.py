import os, torch, argparse, tqdm, random, numpy as np
from util.dataset import create_data
from util.noiser import create_noiser
from util.model import create_model
from torchvision.utils import save_image
from REST_utils import Evaluator, transform_inputs
from eval_utils import get_input_parameters_ibm, read_input_info
from preprocess.grid import Grid
import numpy as np
import pandas as pd
from time import time
from torchvision import transforms
from PIL import Image
from util.RST_prune import Slover
from joblib import Parallel, delayed
from glob import glob

transforms__ = transforms.Compose(
            [
                transforms.Resize(64, interpolation=transforms.InterpolationMode.BILINEAR),
                transforms.CenterCrop(64),
                transforms.Lambda(lambda x: x),
                transforms.ToTensor(),
                transforms.Normalize([0.5], [0.5]),
            ]
        )

device = torch.device("cuda:0")
def seed_everything(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

transforms_ = transforms.Compose(
                [
                    # transforms.ToTensor(),
                    transforms.Normalize([0.5], [0.5]),
                ]
            )
def match_ckpt(ckpt):
    _ckpt = {}
    for k, v in ckpt.items():
        if 'module.' in k:
            k = k.replace('network.module.', 'network.')
        _ckpt[k] = v
    return _ckpt


def route2vertex(route,RGB = True):
    if RGB == True:
        route = route / 2. + 0.5
        route[route <9/10] = 0
        route[route >= 9/10] = 1
        route_up = torch.zeros(route.size())
        route_down = torch.zeros(route.size())
        route_left = torch.zeros(route.size())
        route_right = torch.zeros(route.size())
        route_up[:, :, :, :-1] = route[:,:, :, 1:]
        route_down[:,:, :, 1:] = route[:,:, :, :-1]
        route_left[:,:, :-1, :] = route[:,:, 1:, :]
        route_right[:,:, 1:, :] = route[:,:, :-1, :]

        route = route * torch.sign((route_up + route_down == 1) + \
                                   (route_left + route_right == 1) + \
                                   (route_up + route_down + route_right + route_left == 4))

        route = (route - 0.5) * 2
        return route

def main(args):
    seed_everything(42)
    evaluator = Evaluator(path=os.getcwd() + '/REST_tool/algorithms/libeval.so')
    device = torch.device(f'cuda')

    # prior_set, _, _ = create_data(args.prior, 1, dataset_size=args.num_sample, batch_size=args.batch_size)
    # data_set, _, _ = create_data(args.dataset, 1, dataset_size=args.num_sample, batch_size=args.batch_size)

    noiser = create_noiser(args.noiser, args, device)

    backward_model = create_model(args.method, args, device, noiser, rank=0, direction='b')
    forward_model = create_model(args.method, args, device, noiser, rank=0, direction='f')

    backward_model.to(device)
    forward_model.to(device)

    ckpt = torch.load(args.ckpt, map_location='cpu')
    backward_model.load_state_dict(match_ckpt(ckpt['backward_model']), strict=True)
    forward_model.load_state_dict(match_ckpt(ckpt['forward_model']), strict=True)

    backward_model.eval()
    forward_model.eval()

    if args.benchmark == 'ibm-ada':
        benchmark = ['ibm01.modified.txt', 'ibm02.modified.txt', 'ibm03.modified.txt',
                     'ibm04.modified.txt', 'ibm05.modified.txt', 'adaptec3.dragon70.3d.30.50.90.gr','adaptec4.aplace60.3d.30.50.90.gr']

        batch_size = args.batch_size
        for i in range(len(benchmark)):
            scale = 64
            aligned_generated_routing_twl = 0.
            generated_routing_twl_ = 0.
            generation_time = 0.
            grid_info = read_input_info('preprocess/benchmark/' + benchmark[i])
            grid_param, nets = get_input_parameters_ibm(grid_info)
            env = Grid(grid_param)
            net_array = np.array(
                [[net['net_name'], net['num_pin'], str(net['pins']), net['same_tile']] for net in nets.values()])
            net_pd = pd.DataFrame(net_array, columns=['net_name', 'num_pin', 'pins', 'same_tile'])
            net_pd.num_pin = net_pd.num_pin.astype(int)
            net_pd.same_tile = net_pd.same_tile.apply(lambda x: eval(x))
            net_pd_not_same_tile = net_pd[net_pd.same_tile == False]
            net_pd_not_same_tile = net_pd_not_same_tile.sort_values('num_pin')
            env.planar_cap = torch.tensor(env.planar_cap).to(device)

            all_condition_imgs = []
            all_real_pins = []
            all_random_x = []
            all_random_y = []


            for k in range(len(net_pd_not_same_tile)):
                pins = eval(net_pd_not_same_tile.iloc[k, 2])
                pins = np.array(pins)
                condition_img = torch.zeros((1, 3, scale, scale))

                if max(pins[:, 0]) - min(pins[:, 0]) <= 63 and max(pins[:, 1]) - min(pins[:, 1]) <= 63:
                    random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                              min(pins[:, 0].min(), env.n_x - scale))
                    random_y = random.randint(max(pins[:, 1].max() - scale + 1, 0),
                                              min(pins[:, 1].min(), env.n_y - scale))

                    condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                0]
                    condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                1]
                    for pin in pins:
                        condition_img[0, 0, pin[0] - random_x, pin[1] - random_y] = 255.
                    condition_img /= 255.
                    real_pin = condition_img[:, :1, :, :].detach()
                else:
                    try:
                        random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                                  min(pins[:, 0].min(), env.n_x - scale))
                    except:
                        random_x = min(pins[:, 0].min(), env.n_x - scale)
                    try:
                        random_y = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                                  min(pins[:, 1].min(), env.n_y - scale))
                    except:
                        random_y = min(pins[:, 1].min(), env.n_y - scale)

                    condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                0]
                    condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                1]
                    for pin in pins:
                        try:
                            condition_img[0, 0, pin[0] - random_x, pin[1] - random_y] = 255.
                        except:
                            pass
                    condition_img_copy =  condition_img / 255.
                    real_pin = condition_img_copy[:, :1, :, :].detach()
                    real_pin = real_pin.squeeze()

                all_condition_imgs.append(condition_img)
                all_real_pins.append(real_pin)

                all_random_x.append(random_x)
                all_random_y.append(random_y)
            print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
            with torch.no_grad():
                result_xgens_v1 = [None] * len(all_condition_imgs)  
                start_idx = 0
                while start_idx < len(all_condition_imgs):
                    length_record = []
                    end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                    batch_imgs = []
                    for b_idx in range(start_idx, end_idx):
                        condition_img = transforms_(all_condition_imgs[b_idx]).half()
                        batch_imgs.append(condition_img)
                    batch_imgs = torch.cat(batch_imgs, dim=0).to(device)
                    t0 = time()
                    x_gen_list,_,_,edges_keys  = backward_model.guidance_inferencev1(batch_imgs,args.mode,args.alpha,args.beta,sample=True)
                    t1 = time()
                    x_gen_list = x_gen_list[-1]
                    edges_keys = edges_keys[-1]
                    generation_time += (t1 - t0)
                    for b_idx1 in range(start_idx, end_idx):
                        idx_in_batch = b_idx1 - start_idx
                        result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]
                    for b_id2 in range(x_gen_list.shape[0]):
                        indices = torch.nonzero(batch_imgs[b_id2][0] >= 5 / 6, as_tuple=False)
                        if len(indices) == 0:
                            continue
                        x_gen1 = x_gen_list[b_id2,0,:,:].clone().cpu()
                        x_gen1[x_gen1 <= 9/10] = -1
                        len___ = (x_gen1 > 9/10).sum().item()
                        generated_routing_twl_ += len___
                        if (batch_imgs[b_id2][0] > 5/6).sum().item() == 1:
                            length_record.append(len___)
                            pass
                        else:
                            length_record.append(len___-1)
                            generated_routing_twl_ -= 1
                    x_gen = x_gen_list[:, :1, :, :].detach()
                    x_gen = route2vertex(x_gen)
                    x_vertex = x_gen

                    # t1__ = time()
                    count_ = 0
                    for kk in range(x_vertex.shape[0]):
                        indices = torch.nonzero(batch_imgs[kk][0] >= 5 / 6, as_tuple=False)
                        if len(indices) == 0:
                            continue
                        tc = x_vertex[kk][0].detach().cpu().numpy()
                        tc = np.array(np.where(tc >= 9 / 10)).T
                        # 
                        if(len(tc) >=3000):
                            generated_routing_twl_ -= length_record[count_]
                            continue
                        rx = all_random_x[start_idx + count_]
                        ry = all_random_y[start_idx + count_]
                        tc[:, 0] += rx
                        tc[:, 1] += ry
                        count_ += 1
                        tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                        edges = edges_keys[kk]
                        points = tc
                        point_to_index = {tuple(point): idx for idx, point in enumerate(points)}

                        
                        edge_list_from_edges = []
                        for edge in edges:
                            point1, point2 = edge  
                            u = point_to_index.get(tuple(point1))
                            v = point_to_index.get(tuple(point2))
                            if u is not None and v is not None:
                                edge_list_from_edges.append([u, v])
                            else:
                                pass
                        edge_list = edge_list_from_edges
                        for edge in range(len(edge_list)):
                            u = edge_list[edge][0]
                            v = edge_list[edge][1]
                            if points[u][0] == points[v][0]:
                                env.planar_cap[points[u][0],
                                min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                            elif points[u][1] == points[v][1]:
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]),
                                points[v][1], 0] -= 1
                            else:
                                env.planar_cap[points[u][0],
                                min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]),
                                points[v][1], 0] -= 1

                    print(f'begin:{start_idx},end:{end_idx} generated wl: {aligned_generated_routing_twl}, ' + \
                          f'guidance dsb: {generated_routing_twl_}, ' + \
                          f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                          f'total time: {generation_time}', end='\r')
                    start_idx = end_idx
            print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, ' + \
              f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
              f'total time: {generation_time}')
    if args.benchmark == 'IPSD-s4':
        batch_size = args.batch_size
        generated_routing_twl_ = 0.
        real_routing_twl = 0.
        connected = 0
        generation_time = 0.
        all_condition_imgs = []
        real_routes_len = []
        all_random_x = []
        all_random_y = []
        prior_path = './wlr_a/prior/'
        data_path = './wlr_a/data/'
        c_imgs = sorted(glob(f"{prior_path}/*.png") + glob(f"{prior_path}/*.jpg"))
        r_imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        random_numbers = random.sample(range(75000), 10000)
        for k in random_numbers:
            condition_img = Image.open(c_imgs[k]).convert("RGB")
            real_route = Image.open(r_imgs[k]).convert("RGB")
            real_route = transforms__(real_route)
            route = real_route[0]
            mask__ = route > -0.1
            route_len = mask__.sum().item()
            real_routes_len.append(route_len)
            condition_img = transforms__(condition_img)
            condition_img = condition_img.unsqueeze(0)
            assert condition_img.shape == (1, 3, 64, 64)
            all_condition_imgs.append(condition_img)
            all_random_x.append(0)
            all_random_y.append(0)
        print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
        with torch.no_grad():
            result_xgens_v1 = [None] * len(all_condition_imgs)
            start_idx = 0
            while start_idx < len(all_condition_imgs):
                length_record = []
                end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                batch_imgs = []
                for b_idx in range(start_idx, end_idx):  
                    condition_img = transforms_(all_condition_imgs[b_idx]).half()
                    batch_imgs.append(condition_img)
                batch_imgs = torch.cat(batch_imgs, dim=0).to(device)

                t0 = time()
                x_gen_list, _, _, edges_keys = backward_model.guidance_inferencev1(batch_imgs, sample=True)
                t1 = time()
                x_gen_list = x_gen_list[-1]
                edges_keys = edges_keys[-1]
                # print(f'edges_keys:{edges_keys[0]}')
                generation_time += (t1 - t0)
                for b_idx1 in range(start_idx, end_idx):
                    
                    idx_in_batch = b_idx1 - start_idx
                    result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]
                for b_id2 in range(x_gen_list.shape[0]):
                    real_routing_twl += real_routes_len[start_idx + b_id2]
                    x_gen1 = x_gen_list[b_id2, 0, :, :].clone().cpu()
                    x_gen1[x_gen1 <= 9 / 10] = -1
                    len___ = (x_gen1 > 9 / 10).sum().item()
                    generated_routing_twl_ += len___
                    
                    if (batch_imgs[b_id2][0] > 5 / 6).sum().item() == 1:
                        length_record.append(len___)
                        pass
                    else:
                        length_record.append(len___ - 1)
                        generated_routing_twl_ -= 1
                x_gen = x_gen_list[:, :1, :, :].detach()
                x_gen = route2vertex(x_gen)
                x_vertex = x_gen
                count_ = 0
                for kk in range(x_vertex.shape[0]):
                    tc = x_vertex[kk][0].detach().cpu().numpy()
                    tc = np.array(np.where(tc >= 9 / 10)).T
                    # print(len(tc))
                    if (len(tc) >= 3000):
                        generated_routing_twl_ -= length_record[count_]
                        continue
                    rx = all_random_x[start_idx + count_]
                    ry = all_random_y[start_idx + count_]
                    tc[:, 0] += rx
                    tc[:, 1] += ry
                    count_ += 1
                    tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                    edges = edges_keys[kk]
                    points = tc
                    point_to_index = {tuple(point): idx for idx, point in enumerate(points)}
                    edge_list_from_edges = []
                    for edge in edges:
                        point1, point2 = edge
                        u = point_to_index.get(tuple(point1))
                        v = point_to_index.get(tuple(point2))
                        if u is not None and v is not None:
                            edge_list_from_edges.append([u, v])
                        else:
                            pass

                print(f'begin:{start_idx},end:{end_idx}, ' + \
                      f'guidance routing wl: {generated_routing_twl_}, real routing wl:{real_routing_twl}, ' + \
                      f'ratio: {generated_routing_twl_ / real_routing_twl}, ' + \
                      # f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                      f'total time: {generation_time}', end='\r')
                start_idx = end_idx

            print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
              f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
              f'total time: {generation_time}')
    if args.benchmark == 'IPSD-s':
        batch_size = args.batch_size
        generated_routing_twl_ = 0.
        real_routing_twl = 0.
        connected = 0
        generation_time = 0.
        all_condition_imgs = []
        real_routes_len = []
        all_random_x = []
        all_random_y = []
        prior_path = './wlr_b/prior/'
        data_path = './wlr_b/data/'
        c_imgs = sorted(glob(f"{prior_path}/*.png") + glob(f"{prior_path}/*.jpg"))
        r_imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        random_numbers = random.sample(range(75000), 10000)
        for k in random_numbers:
            condition_img = Image.open(c_imgs[k]).convert("RGB")
            real_route = Image.open(r_imgs[k]).convert("RGB")
            real_route = transforms__(real_route)
            route = real_route[0]
            mask__ = route > -0.1
            route_len = mask__.sum().item()
            real_routes_len.append(route_len)
            condition_img = transforms__(condition_img)
            condition_img = condition_img.unsqueeze(0)
            assert condition_img.shape == (1, 3, 64, 64)
   
            all_condition_imgs.append(condition_img)
           
            all_random_x.append(0)
            all_random_y.append(0)
        print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
        with torch.no_grad():
            result_xgens_v1 = [None] * len(all_condition_imgs)  
            start_idx = 0
            while start_idx < len(all_condition_imgs):
                length_record = []
                end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                batch_imgs = []
                for b_idx in range(start_idx, end_idx):
                    condition_img = transforms_(all_condition_imgs[b_idx]).half()
                    batch_imgs.append(condition_img)
                batch_imgs = torch.cat(batch_imgs, dim=0).to(device)

                t0 = time()
                x_gen_list, _, _, edges_keys = backward_model.guidance_inferencev1(batch_imgs, sample=True)
                t1 = time()
                x_gen_list = x_gen_list[-1]
                edges_keys = edges_keys[-1]
                generation_time += (t1 - t0)
                for b_idx1 in range(start_idx, end_idx):
                    idx_in_batch = b_idx1 - start_idx
                    result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]
                for b_id2 in range(x_gen_list.shape[0]):
                    real_routing_twl += real_routes_len[start_idx + b_id2]
                    x_gen1 = x_gen_list[b_id2, 0, :, :].clone().cpu()
                    x_gen1[x_gen1 <= 9 / 10] = -1
                    len___ = (x_gen1 > 9 / 10).sum().item()
                    generated_routing_twl_ += len___
                    if (batch_imgs[b_id2][0] > 5 / 6).sum().item() == 1:
                        length_record.append(len___)
                        pass
                    else:
                        length_record.append(len___ - 1)
                        generated_routing_twl_ -= 1
                x_gen = x_gen_list[:, :1, :, :].detach()
                x_gen = route2vertex(x_gen)
                x_vertex = x_gen
                count_ = 0
                for kk in range(x_vertex.shape[0]):
                    tc = x_vertex[kk][0].detach().cpu().numpy()
                    tc = np.array(np.where(tc >= 9 / 10)).T
                    if (len(tc) >= 3000):
                        generated_routing_twl_ -= length_record[count_]
                        continue
                    rx = all_random_x[start_idx + count_]
                    ry = all_random_y[start_idx + count_]
                    tc[:, 0] += rx
                    tc[:, 1] += ry
                    count_ += 1
                    tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                    edges = edges_keys[kk]
                    points = tc
                    point_to_index = {tuple(point): idx for idx, point in enumerate(points)}
                    edge_list_from_edges = []
                    for edge in edges:
                        point1, point2 = edge 
                        u = point_to_index.get(tuple(point1))
                        v = point_to_index.get(tuple(point2))
                        if u is not None and v is not None:
                            edge_list_from_edges.append([u, v])
                        else:
                            pass

                print(f'begin:{start_idx},end:{end_idx}, ' + \
                      f'guidance routing wl: {generated_routing_twl_}, real routing wl:{real_routing_twl}, ' + \
                      f'ratio: {generated_routing_twl_ / real_routing_twl}, ' + \
                      # f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                      f'total time: {generation_time}', end='\r')
                start_idx = end_idx
            print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
              f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
              f'total time: {generation_time}')
    if args.benchmark == 'IPSD-b4':
        batch_size = args.batch_size
        generated_routing_twl_ = 0.
        real_routing_twl = 0.
        connected = 0
        generation_time = 0.
        postprocess_time = 0.
        all_condition_imgs = []
        real_routes_len = []
        all_random_x = []
        all_random_y = []
        prior_path = './wlr_c/prior/'
        data_path = './wlr_c/data/'
        c_imgs = sorted(glob(f"{prior_path}/*.png") + glob(f"{prior_path}/*.jpg"))
        r_imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        random_numbers = random.sample(range(75000), 10000)
        for k in random_numbers:
            condition_img = Image.open(c_imgs[k]).convert("RGB")
            real_route = Image.open(r_imgs[k]).convert("RGB")
            real_route = transforms__(real_route)
            route = real_route[0]
            mask__ = route > -0.1
            route_len = mask__.sum().item()
            real_routes_len.append(route_len)
            condition_img = transforms__(condition_img)
            condition_img = condition_img.unsqueeze(0)
            assert condition_img.shape == (1, 3, 64, 64)
            all_condition_imgs.append(condition_img)
            all_random_x.append(0)
            all_random_y.append(0)
        print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
        with torch.no_grad():
            result_xgens_v1 = [None] * len(all_condition_imgs)
            start_idx = 0
            while start_idx < len(all_condition_imgs):
                length_record = []
                end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                batch_imgs = []
                for b_idx in range(start_idx, end_idx):
                    condition_img = transforms_(all_condition_imgs[b_idx]).half()
                    batch_imgs.append(condition_img)
                batch_imgs = torch.cat(batch_imgs, dim=0).to(device)

                t0 = time()
                x_gen_list, _, _, edges_keys = backward_model.guidance_inferencev1(batch_imgs, sample=True)
                t1 = time()
                x_gen_list = x_gen_list[-1]
                edges_keys = edges_keys[-1]
                generation_time += (t1 - t0)
                for b_idx1 in range(start_idx, end_idx):
                    idx_in_batch = b_idx1 - start_idx
                    result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]
                for b_id2 in range(x_gen_list.shape[0]):
                    real_routing_twl += real_routes_len[start_idx + b_id2]
                    x_gen1 = x_gen_list[b_id2, 0, :, :].clone().cpu()
                    x_gen1[x_gen1 <= 9 / 10] = -1
                    len___ = (x_gen1 > 9 / 10).sum().item()
                    generated_routing_twl_ += len___
                    if (batch_imgs[b_id2][0] > 5 / 6).sum().item() == 1:
                        length_record.append(len___)
                        pass
                    else:
                        length_record.append(len___ - 1)
                        generated_routing_twl_ -= 1
                x_gen = x_gen_list[:, :1, :, :].detach()
                x_gen = route2vertex(x_gen)
                x_vertex = x_gen
                count_ = 0
                for kk in range(x_vertex.shape[0]):
                    tc = x_vertex[kk][0].detach().cpu().numpy()
                    tc = np.array(np.where(tc >= 9 / 10)).T
                    if (len(tc) >= 3000):
                        generated_routing_twl_ -= length_record[count_]
                        continue
                    rx = all_random_x[start_idx + count_]
                    ry = all_random_y[start_idx + count_]
                    tc[:, 0] += rx
                    tc[:, 1] += ry
                    count_ += 1
                    tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                    edges = edges_keys[kk]
                    points = tc
                    point_to_index = {tuple(point): idx for idx, point in enumerate(points)}

                    edge_list_from_edges = []
                    for edge in edges:
                        point1, point2 = edge
                        u = point_to_index.get(tuple(point1))
                        v = point_to_index.get(tuple(point2))
                        if u is not None and v is not None:
                            edge_list_from_edges.append([u, v])
                        else:
                            pass

                print(f'begin:{start_idx},end:{end_idx}, ' + \
                      f'guidance routing wl: {generated_routing_twl_}, real routing wl:{real_routing_twl}, ' + \
                      f'ratio: {generated_routing_twl_ / real_routing_twl}, ' + \
                      # f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                      f'total time: {generation_time}', end='\r')
                start_idx = end_idx
            print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
              f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
              f'total time: {generation_time}')
    if args.benchmark == 'IPSD-b':
        batch_size = args.batch_size
        generated_routing_twl_ = 0.
        real_routing_twl = 0.
        connected = 0
        generation_time = 0.
        postprocess_time = 0.
        all_condition_imgs = []
        real_routes_len = []
        all_random_x = []
        all_random_y = []
        prior_path = './wlr_d/prior/'
        data_path = './wlr_d/data/'
        c_imgs = sorted(glob(f"{prior_path}/*.png") + glob(f"{prior_path}/*.jpg"))
        r_imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        random_numbers = random.sample(range(75000), 10000)
        for k in random_numbers:
            condition_img = Image.open(c_imgs[k]).convert("RGB")
            real_route = Image.open(r_imgs[k]).convert("RGB")
            real_route = transforms__(real_route)
            route = real_route[0]
            mask__ = route > -0.1
            route_len = mask__.sum().item()
            real_routes_len.append(route_len)
            condition_img = transforms__(condition_img)
            condition_img = condition_img.unsqueeze(0)
            assert condition_img.shape == (1, 3, 64, 64)
            all_condition_imgs.append(condition_img)
            all_random_x.append(0)
            all_random_y.append(0)
        print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
        with torch.no_grad():
            result_xgens_v1 = [None] * len(all_condition_imgs)
            start_idx = 0
            while start_idx < len(all_condition_imgs):
                length_record = []
                end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                batch_imgs = []
                for b_idx in range(start_idx, end_idx):
                    condition_img = transforms_(all_condition_imgs[b_idx]).half()
                    batch_imgs.append(condition_img)
                batch_imgs = torch.cat(batch_imgs, dim=0).to(device)

                t0 = time()
                x_gen_list, _, _, edges_keys = backward_model.guidance_inferencev1(batch_imgs, sample=True)
                t1 = time()
                x_gen_list = x_gen_list[-1]
                edges_keys = edges_keys[-1]
                generation_time += (t1 - t0)
                for b_idx1 in range(start_idx, end_idx):
                    idx_in_batch = b_idx1 - start_idx
                    result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]
                for b_id2 in range(x_gen_list.shape[0]):
                    real_routing_twl += real_routes_len[start_idx + b_id2]
                    x_gen1 = x_gen_list[b_id2, 0, :, :].clone().cpu()
                    x_gen1[x_gen1 <= 9 / 10] = -1
                    len___ = (x_gen1 > 9 / 10).sum().item()
                    generated_routing_twl_ += len___
                    if (batch_imgs[b_id2][0] > 5 / 6).sum().item() == 1:
                        length_record.append(len___)
                        pass
                    else:
                        length_record.append(len___ - 1)
                        generated_routing_twl_ -= 1
                x_gen = x_gen_list[:, :1, :, :].detach()
                x_gen = route2vertex(x_gen)
                x_vertex = x_gen
                count_ = 0
                for kk in range(x_vertex.shape[0]):
                    tc = x_vertex[kk][0].detach().cpu().numpy()
                    tc = np.array(np.where(tc >= 9 / 10)).T
                    if (len(tc) >= 3000):
                        generated_routing_twl_ -= length_record[count_]
                        continue
                    rx = all_random_x[start_idx + count_]
                    ry = all_random_y[start_idx + count_]
                    tc[:, 0] += rx
                    tc[:, 1] += ry
                    count_ += 1
                    tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                    edges = edges_keys[kk]
                    points = tc
                    point_to_index = {tuple(point): idx for idx, point in enumerate(points)}
                    edge_list_from_edges = []
                    for edge in edges:
                        point1, point2 = edge
                        u = point_to_index.get(tuple(point1))
                        v = point_to_index.get(tuple(point2))
                        if u is not None and v is not None:
                            edge_list_from_edges.append([u, v])
                        else:
                            pass

                print(f'begin:{start_idx},end:{end_idx}, ' + \
                      f'guidance routing wl: {generated_routing_twl_}, real routing wl:{real_routing_twl}, ' + \
                      f'ratio: {generated_routing_twl_ / real_routing_twl}, ' + \
                      # f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                      f'total time: {generation_time}', end='\r')
                start_idx = end_idx
                # if start_idx >=63:
                #     break
            print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
              f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
              f'total time: {generation_time}')

    if args.benchmark == 'abl-NN':
        benchmark = ['ibm01.modified.txt']
        batch_size = args.batch_size
        for i in range(len(benchmark)):
            scale = 64
            aligned_generated_routing_twl = 0.
            generated_routing_twl_ = 0.
            generation_time = 0.
            connected = 0
            grid_info = read_input_info('preprocess/benchmark/' + benchmark[i])
            grid_param, nets = get_input_parameters_ibm(grid_info)
            env = Grid(grid_param)
            net_array = np.array(
                [[net['net_name'], net['num_pin'], str(net['pins']), net['same_tile']] for net in nets.values()])
            net_pd = pd.DataFrame(net_array, columns=['net_name', 'num_pin', 'pins', 'same_tile'])
            net_pd.num_pin = net_pd.num_pin.astype(int)
            net_pd.same_tile = net_pd.same_tile.apply(lambda x: eval(x))
            net_pd_not_same_tile = net_pd[net_pd.same_tile == False]
            net_pd_not_same_tile = net_pd_not_same_tile.sort_values('num_pin')
            env.planar_cap = torch.tensor(env.planar_cap).to(device)

            all_condition_imgs = []
            all_real_pins = []
            all_pins_data = []
            all_random_x = []
            all_random_y = []


            for k in range(len(net_pd_not_same_tile)):
                pins = eval(net_pd_not_same_tile.iloc[k, 2])
                pins = np.array(pins)
                condition_img = torch.zeros((1, 3, scale, scale))

                if max(pins[:, 0]) - min(pins[:, 0]) <= 63 and max(pins[:, 1]) - min(pins[:, 1]) <= 63:
                    random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                              min(pins[:, 0].min(), env.n_x - scale))
                    random_y = random.randint(max(pins[:, 1].max() - scale + 1, 0),
                                              min(pins[:, 1].min(), env.n_y - scale))

                    condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                0]
                    condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                1]
                    for pin in pins:
                        condition_img[0, 0, pin[0] - random_x, pin[1] - random_y] = 255.
                    condition_img /= 255.
                    real_pin = condition_img[:, :1, :, :].detach()
                else:
                    try:
                        random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                                  min(pins[:, 0].min(), env.n_x - scale))
                    except:
                        random_x = min(pins[:, 0].min(), env.n_x - scale)
                    try:
                        random_y = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                                  min(pins[:, 1].min(), env.n_y - scale))
                    except:
                        random_y = min(pins[:, 1].min(), env.n_y - scale)

                    condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                0]
                    condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                1]
                    for pin in pins:
                        try:
                            condition_img[0, 0, pin[0] - random_x, pin[1] - random_y] = 255.
                        except:
                            pass
                    condition_img_copy =  condition_img / 255.
                    real_pin = condition_img_copy[:, :1, :, :].detach()
                    real_pin = real_pin.squeeze()
                all_condition_imgs.append(condition_img)
                all_real_pins.append(real_pin)
                # all_pins_data.append(pins)
                all_random_x.append(random_x)
                all_random_y.append(random_y)
            print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
            with torch.no_grad():
                result_xgens_v1 = [None] * len(all_condition_imgs) 
                result_need_second_pass = []  
                start_idx = 0
                while start_idx < len(all_condition_imgs):
                    length_record = []
                    end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                    batch_imgs = []
                    for b_idx in range(start_idx, end_idx):
                        condition_img = transforms_(all_condition_imgs[b_idx]).half()
                        batch_imgs.append(condition_img)
                    batch_imgs = torch.cat(batch_imgs, dim=0).to(device)

                    t0 = time()
                    slover = Slover(args={'batch_size': batch_imgs.shape[0], 'mode': 'of'})
                    bg = batch_imgs.clone()
                    x_other = batch_imgs.clone()
                    c_raw = batch_imgs[:, 0, :, :]
                    c_raw = (c_raw >= 5 / 6).to(torch.float32)
                    slover.SetTarget(bg, pins=x_other, conditions=c_raw, steinerpoints=None)
                    results = Parallel(n_jobs=slover.args['batch_size'])(
                        delayed(slover.Parallelization_tasks)(i) for i in range(slover.args['batch_size'])
                    )
                    mask_route = torch.full((x_other.shape[0], 3, 64, 64), -1.0, dtype=torch.float32).to("cpu")
                    mask_route[:, 1:, :, :] = bg[:, 1:, :, :]
                    mask_route[:, 2:, :, :] = bg[:, 2:, :, :]

                    keys_ = {}
                    for idx, res in enumerate(results):
                        mask_route[idx][0], keys = res[0][0], res[1]
                        keys_[idx] = keys
                    x_gen_list = mask_route
                    t1 = time()
                    x_gen_list = mask_route
                    edges_keys = keys_
                    generation_time += (t1 - t0)

                    for b_idx1 in range(start_idx, end_idx):
                        idx_in_batch = b_idx1 - start_idx
                        result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]

                    for b_id2 in range(x_gen_list.shape[0]):
                        x_gen1 = x_gen_list[b_id2,0,:,:].clone().cpu()
                        x_gen1[x_gen1 <= 9/10] = -1
                        len___ = (x_gen1 > 9/10).sum().item()
                        generated_routing_twl_ += len___
                        if (batch_imgs[b_id2][0] > 5/6).sum().item() == 1:
                            length_record.append(len___)
                            pass
                        else:
                            length_record.append(len___-1)
                            generated_routing_twl_ -= 1
                    x_gen = x_gen_list[:, :1, :, :].detach()
                    x_gen = route2vertex(x_gen)

                    x_vertex = x_gen


                    count_ = 0
                    for kk in range(x_vertex.shape[0]):
                        tc = x_vertex[kk][0].detach().cpu().numpy()
                        tc = np.array(np.where(tc >= 9 / 10)).T
                        if(len(tc) >=3000):
                            generated_routing_twl_ -= length_record[count_]
                            continue
                        rx = all_random_x[start_idx + count_]
                        ry = all_random_y[start_idx + count_]
                        tc[:, 0] += rx
                        tc[:, 1] += ry
                        count_ += 1
                        tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                        edges = edges_keys[kk]
                        points = tc
                        point_to_index = {tuple(point): idx for idx, point in enumerate(points)}
                        edge_list_from_edges = []
                        for edge in edges:
                            point1, point2 = edge
                            u = point_to_index.get(tuple(point1))
                            v = point_to_index.get(tuple(point2))
                            if u is not None and v is not None:
                                edge_list_from_edges.append([u, v])
                            else:
                                pass
                        edge_list = edge_list_from_edges
                        for edge in range(len(edge_list)):
                            u = edge_list[edge][0]
                            v = edge_list[edge][1]
                            if points[u][0] == points[v][0]:
                                env.planar_cap[points[u][0],
                                min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                            elif points[u][1] == points[v][1]:
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]),
                                points[v][1], 0] -= 1
                            else:
                                env.planar_cap[points[u][0],
                                min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]),
                                points[v][1], 0] -= 1

                    print(f'begin:{start_idx},end:{end_idx} generated wl: {aligned_generated_routing_twl}, ' + \
                          f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
                          f'generation time: {generation_time}s, postprocess_time: {postprocess_time}s, ' + \
                          f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                          f'total time: {generation_time + postprocess_time}', end='\r')
                    start_idx = end_idx
            print()

        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
                  f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                  f'total time: {generation_time}')
    if args.benchmark == 'abl-RST':
        benchmark = ['ibm01.modified.txt']
        batch_size = args.batch_size
        for i in range(len(benchmark)):
            scale = 64
            aligned_generated_routing_twl = 0.
            generated_routing_twl_ = 0.
            generation_time = 0.
            connected = 0
            grid_info = read_input_info('preprocess/benchmark/' + benchmark[i])
            grid_param, nets = get_input_parameters_ibm(grid_info)
            env = Grid(grid_param)
            net_array = np.array(
                [[net['net_name'], net['num_pin'], str(net['pins']), net['same_tile']] for net in nets.values()])
            net_pd = pd.DataFrame(net_array, columns=['net_name', 'num_pin', 'pins', 'same_tile'])
            net_pd.num_pin = net_pd.num_pin.astype(int)
            net_pd.same_tile = net_pd.same_tile.apply(lambda x: eval(x))
            # net_pd_same_tile = net_pd[net_pd.same_tile == True]
            net_pd_not_same_tile = net_pd[net_pd.same_tile == False]
            net_pd_not_same_tile = net_pd_not_same_tile.sort_values('num_pin')
            env.planar_cap = torch.tensor(env.planar_cap).to(device)

            all_condition_imgs = []
            all_real_pins = []
            all_random_x = []
            all_random_y = []


            for k in range(len(net_pd_not_same_tile)):
                pins = eval(net_pd_not_same_tile.iloc[k, 2])
                pins = np.array(pins)
                condition_img = torch.zeros((1, 3, scale, scale))

                if max(pins[:, 0]) - min(pins[:, 0]) <= 63 and max(pins[:, 1]) - min(pins[:, 1]) <= 63:
                    random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                              min(pins[:, 0].min(), env.n_x - scale))
                    random_y = random.randint(max(pins[:, 1].max() - scale + 1, 0),
                                              min(pins[:, 1].min(), env.n_y - scale))

                    condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                0]
                    condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                1]
                    for pin in pins:
                        condition_img[0, 0, pin[0] - random_x, pin[1] - random_y] = 255.
                    condition_img /= 255.
                    real_pin = condition_img[:, :1, :, :].detach()
                else:
                    try:
                        random_x = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                                  min(pins[:, 0].min(), env.n_x - scale))
                    except:
                        random_x = min(pins[:, 0].min(), env.n_x - scale)
                    try:
                        random_y = random.randint(max(pins[:, 0].max() - scale + 1, 0),
                                                  min(pins[:, 1].min(), env.n_y - scale))
                    except:
                        random_y = min(pins[:, 1].min(), env.n_y - scale)

                    condition_img[0, 1, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                0]
                    condition_img[0, 2, :, :] = env.planar_cap[random_x:random_x + scale, random_y:random_y + scale,
                                                1]
                    for pin in pins:
                        try:
                            condition_img[0, 0, pin[0] - random_x, pin[1] - random_y] = 255.
                        except:
                            pass
                    condition_img_copy =  condition_img / 255.
                    real_pin = condition_img_copy[:, :1, :, :].detach()
                    real_pin = real_pin.squeeze()
                all_condition_imgs.append(condition_img)
                all_real_pins.append(real_pin)
                # all_pins_data.append(pins)
                all_random_x.append(random_x)
                all_random_y.append(random_y)
            print(f'len of net_pd_not_same_tile: {len(all_condition_imgs)}')
            with torch.no_grad():
                result_xgens_v1 = [None] * len(all_condition_imgs)
                result_need_second_pass = []
                start_idx = 9128
                while start_idx < len(all_condition_imgs):
                    length_record = []
                    end_idx = min(start_idx + batch_size, len(all_condition_imgs))
                    batch_imgs = []
                    for b_idx in range(start_idx, end_idx):
                        condition_img = transforms_(all_condition_imgs[b_idx]).half()
                        batch_imgs.append(condition_img)
                    batch_imgs = torch.cat(batch_imgs, dim=0).to(device)

                    t0 = time()
                    x_gen_list,_,_  = backward_model.guidance_inference_abl(batch_imgs,sample=True)
                    t1 = time()
                    x_gen_list = x_gen_list[-1]
                    generation_time += (t1 - t0)

                    for b_idx1 in range(start_idx, end_idx):
                        idx_in_batch = b_idx1 - start_idx
                        result_xgens_v1[b_idx1] = x_gen_list[idx_in_batch][-1]

                    for b_id2 in range(x_gen_list.shape[0]):
                        x_gen1 = x_gen_list[b_id2,0,:,:].clone().cpu()
                        x_gen1[x_gen1 <= 9/10] = -1
                        len___ = (x_gen1 > 9/10).sum().item()
                        generated_routing_twl_ += len___
                        if (batch_imgs[b_id2][0] > 5/6).sum().item() == 1:
                            length_record.append(len___)
                            pass
                        else:
                            length_record.append(len___-1)
                            generated_routing_twl_ -= 1
                    x_gen = x_gen_list[:, :1, :, :].detach()
                    x_gen = route2vertex(x_gen)

                    real_pin_ = torch.cat(all_real_pins[start_idx:end_idx], dim=0)
                    for b_idx3 in range(start_idx, end_idx):

                        pin_index = torch.where(all_real_pins[b_idx3] == 1.)
                        mask_route = torch.zeros_like(all_real_pins[b_idx3])
                        mask_route[0, pin_index[1].min(): pin_index[1].max() + 1,
                        pin_index[2].min(): pin_index[2].max() + 1] = 1.
                        x_gen[b_idx3-start_idx] = x_gen[b_idx3-start_idx] * mask_route
                    x_vertex = torch.max(x_gen, real_pin_)

                    count_ = 0
                    for kk in range(x_vertex.shape[0]):
                        tc = x_vertex[kk][0].detach().cpu().numpy()
                        tc = np.array(np.where(tc >= 5 / 6)).T
        
                        rx = all_random_x[start_idx + count_]
                        ry = all_random_y[start_idx + count_]
                        tc[:, 0] += rx
                        tc[:, 1] += ry
                        count_ += 1
                        tc = np.array(list(set([tuple(cnnt) for cnnt in tc])))
                        gst_length, sp_list, edge_list = evaluator.gst_rsmt(tc)
                        aligned_generated_routing_twl += gst_length

                        sp_list = np.int32(sp_list)
                        if len(sp_list) != 0:
                            points = np.concatenate([tc, sp_list], 0)
                        else:
                            points = tc
                        for edge in range(len(edge_list)):
                            u = edge_list[edge][0]
                            v = edge_list[edge][1]
                            if points[u][0] == points[v][0]:
                                env.planar_cap[points[u][0],
                                min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                            elif points[u][1] == points[v][1]:
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]),
                                points[v][1], 0] -= 1
                            else:
                                env.planar_cap[points[u][0],
                                min(points[u][1], points[v][1]): max(points[u][1], points[v][1]), 1] -= 1
                                env.planar_cap[min(points[u][0], points[v][0]): max(points[u][0], points[v][0]),
                                points[v][1], 0] -= 1

                    print(f'begin:{start_idx},end:{end_idx} generated wl: {aligned_generated_routing_twl}, ' + \
                          f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
                          f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                          f'total time: {generation_time}', end='\r')
                    start_idx = end_idx
            print()
        print(f'net: {k}/{len(net_pd_not_same_tile)}, generated wl: {aligned_generated_routing_twl}, ' + \
              f'guidance dsb: {generated_routing_twl_}, connected: {connected}, ' + \
                  f'overflow: {(-env.planar_cap * (env.planar_cap < 0)).sum().item()} ' + \
                  f'total time: {generation_time}')

def create_parser():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--num_sample', type=int, default=128, help='number of samples')
    argparser.add_argument('--batch_size', type=int, default=64, help='batch size')

    argparser.add_argument('--method', type=str, default='dsb', help='method')
    argparser.add_argument('--simplify', action='store_true', help='whether to use simplified DSB')
    argparser.add_argument('--reparam', type=str, default='term',
                           help='whether to use reparameterized DSB, "term" for TR-DSB, "flow" for FR-DSB')
    argparser.add_argument('--noiser', type=str, default='flow',
                           help='noiser type, "flow" noiser for Flow Matching models, "dsb" noiser for DSB models')
    argparser.add_argument('--gamma_type', type=str, default='constant', help='gamma schedule for DSB')
    argparser.add_argument('--training_timesteps', type=int, default=24, help='training timesteps')
    argparser.add_argument('--inference_timesteps', type=int, default=24, help='inference timesteps')

    argparser.add_argument('--network', type=str, default='uvit-b', help='network architecture to use')
    argparser.add_argument('--use_amp', action='store_true', help='whether to use mixed-precision training')

    # argparser.add_argument('--prior', type=str, default='standard', help='prior distribution')
    # argparser.add_argument('--dataset', type=str, default='checkerboard:4', help='data distribution')

    argparser.add_argument('--exp_name', type=str, default='nthu_inference', help='name of experiment')
    argparser.add_argument('--ckpt', type=str, default=None, help='checkpoint to load')

    argparser.add_argument('--benchmark', type=str, choices=['ibm-ada','abl-NN','abl-RST','IPSD-s4','IPSD-s','IPSD-b4','IPSD-b'], default='ibm-ada', help='choose benchmark')
    argparser.add_argument('--scale', type=str, default='64', help='scale of the benchmark')

    
    argparser.add_argument('--mode', type=str, choices=['of','wl','avg_of','mix'], default='of', help='choose the gradient search mode')
   
    argparser.add_argument('--alpha', type=float, default=None, help='choose tune alpha')
    argparser.add_argument('--beta', type=float, default=None, help='choose tune beta')

    return argparser


if __name__ == '__main__':
    # with concurrent.futures.ProcessPoolExecutor(max_workers=slover.args['batch_size']) as executor:

    argparser = create_parser()
    args = argparser.parse_args()
    # executor = concurrent.futures.ProcessPoolExecutor(max_workers=args.batch_size)

    if 'dsb' in args.method:
        assert args.training_timesteps == args.inference_timesteps

    main(args)
