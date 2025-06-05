import os.path
import torch
from functools import partial
from util.network import ResMLP, DhariwalUNet
from util.vision_transformer import VisionTransformer
from util.RST_prune import *
import concurrent.futures
from joblib import Parallel, delayed
from torchvision import transforms

def check_size(x, coef):
    if isinstance(coef, (int, float)):
        return coef
    elif isinstance(coef, dict):
        for k, v in coef.items():
            if isinstance(v, torch.Tensor):
                while len(v.shape) < len(x.shape):
                    v = v.unsqueeze(-1)
                coef[k] = v
    elif isinstance(coef, torch.Tensor):
        while len(coef.shape) < len(x.shape):
            coef = coef.unsqueeze(-1)
    return coef


class BaseModel(torch.nn.Module):
    def __init__(self, args, device, noiser, rank):
        super().__init__()
        self.args = args
        self.device = device
        self.noiser = noiser
        self.rank = rank

        if self.args.network == 'mlp':
            self.network = ResMLP(dim_in=2, dim_out=2, dim_hidden=128, num_layers=5, n_cond=self.noiser.training_timesteps)
        else: 
            img_resolution = 32
            if "-512" in args.dataset:
                img_resolution = 64
            if "eda" in args.dataset:
                img_resolution = 64    
            in_channels = 3
            if self.args.dataset == "afhq-cat-64" or self.args.dataset == "celeba-64":
                in_channels = 3

            if args.network == 'adm':
                self.network = DhariwalUNet(
                    img_resolution=img_resolution,
                    in_channels=in_channels,
                    out_channels=in_channels,
                    label_dim=0,
                    augment_dim=0,
                    model_channels=96,
                    channel_mult=[1,2,2,2],
                    channel_mult_emb=4,
                    num_blocks=4,
                    attn_resolutions=[8, 4],
                    dropout=0.10,
                    label_dropout=0,
                )
            elif args.network == 'uvit-b':
                self.network = VisionTransformer(
                    img_size=64,
                    patch_size=4, embed_dim=512, 
                    num_heads=8, mlp_ratio=4, qkv_bias=True,
                    depth=13,
                    norm_layer=partial(torch.nn.LayerNorm, eps=1e-6), 
                    num_classes=-1,
                    use_fp16=self.args.use_amp,
                )

        self.network.to(self.device)
        self.network = torch.nn.parallel.DistributedDataParallel(self.network, device_ids=[self.device], output_device=self.device) if (hasattr(self.args, 'gpus') and self.args.gpus > 1) else self.network

    def target(self, x_0, x_1, x_t, t):
        raise NotImplementedError
    def forward(self, x_t, t):
        t = self.noiser.timestep_map[t]
        x = self.network(x_t, t)
        return x

    def predict_boundary(self, x_0, x_t, t):
        raise NotImplementedError
    
    def predict_next(self, x_0, x_t, t):
        x_0, x_1 = self.predict_boundary(x_0, x_t, t)
        x_t = self.noiser(x_0, x_1, t + 1)
        return x_t

    def inference(self, x_0, return_all=False):
        self.eval()
        x_t_all = [x_0.clone()]
        with torch.no_grad():
            x_t = x_0
            ones = torch.ones(size=(x_t.shape[0],), dtype=torch.int64, device=x_t.device)
            for t in range(self.noiser.num_timesteps):
                with torch.autocast(device_type="cuda", enabled=self.args.use_amp):
                    x_t = self.predict_next(x_0, x_t, ones * t)
                x_t = x_t.float()
                if return_all:
                    x_t_all.append(x_t.clone())
        if return_all:
            return x_t, torch.stack(x_t_all, dim=0)
        else:
            return x_t


class Diffusion(BaseModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.prediction_type = self.args.method.split(':')[1] if ':' in self.args.method else 'x1'

    def target(self, x_0, x_1, x_t, t):
        if self.prediction_type == 'x1':
            return x_1
        elif self.prediction_type == 'x0':
            return x_0
        elif self.prediction_type == 'v':
            coef = check_size(x_t, self.noiser.coefficient(t))
            v = coef['coef1'] * x_0 - coef['coef0'] * x_1
            return v

    def predict_boundary(self, x_0, x_t, t):
        coef_t = check_size(x_t, self.noiser.coefficient(t))
        if self.prediction_type == 'x1':
            x_1 = self.forward(x_t, t)
            x_0 = (x_t - coef_t['coef1'] * x_1) / coef_t['coef0']
        elif self.prediction_type == 'x0':
            x_0 = self.forward(x_t, t)
            x_1 = (x_t - coef_t['coef0'] * x_0) / coef_t['coef1']
        elif self.prediction_type == 'v':
            v = self.forward(x_t, t)
            x_0 = coef_t['coef0'] * x_t + coef_t['coef1'] * v
            x_1 = coef_t['coef1'] * x_t - coef_t['coef0'] * v
        return x_0, x_1


class Flow(BaseModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def target(self, x_0, x_1, x_t, t):
        return x_1 - x_0
    
    def predict_next(self, x_0, x_t, t):
        coef_t_1 = check_size(x_t, self.noiser.coefficient(t))['coef1']
        coef_t_plus_one_1 = check_size(x_t, self.noiser.coefficient(t + 1))['coef1']
        v_pred = self.forward(x_t, t)
        x_t = x_t + (coef_t_plus_one_1 - coef_t_1) * v_pred
        return x_t


class DSB(BaseModel):
    def __init__(self, *args, direction=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.direction = direction
        self.num_timesteps = self.noiser.training_timesteps

        self.noiser.prepare_gamma_dsb()
        self.gammas = self.noiser.gammas

    def get_coef_ts(self, x, t, delta=1):
        coef_t = check_size(x, self.noiser.coefficient(t))
        coef_t_other = check_size(x, self.noiser.coefficient(t + delta))
        return coef_t, coef_t_other

    def _forward(self, x, t):
        x_other = self.forward(x, t)
        if self.args.reparam == 'flow':
            v_pred = x_other
            if self.direction == 'b':
                coef_t, coef_t_next = self.get_coef_ts(x, t, 1)
                x = x + (coef_t_next['coef1'] - coef_t['coef1']) * v_pred
            elif self.direction == 'f':
                coef_t, coef_t_next = self.get_coef_ts(x, self.num_timesteps - t, -1)
                x = x + (coef_t_next['coef0'] - coef_t['coef0']) * v_pred
        elif self.args.reparam == 'term':
            if self.direction == 'b':
                coef_t, coef_t_next = self.get_coef_ts(x, t, 1)
                x_1 = x_other
                x_0 = (x - coef_t['coef1'] * x_1) / coef_t['coef0']
            elif self.direction == 'f':
                coef_t, coef_t_next = self.get_coef_ts(x, self.num_timesteps - t, -1)
                x_0 = x_other
                x_1 = (x - coef_t['coef0'] * x_0) / coef_t['coef1']
            x = coef_t_next['coef0'] * x_0 + coef_t_next['coef1'] * x_1
        else:
            x = x_other
        return x


    def guidance_fowardv1(self,x,t,c_raw,bg,count_,slover: Slover):
        keys_ = {}
        x_other = self.forward(x, t)
        points_ = (c_raw >= 5 / 6).nonzero(as_tuple=True)
        x_other[:, 0, :, :][points_] = 1

        c_raw_ = c_raw.to("cpu")
        slover.SetTarget( bg,pins=x_other, conditions=c_raw_,steinerpoints=None)

        results = Parallel(n_jobs=slover.args['batch_size'])(
            delayed(slover.Parallelization_tasks)(i) for i in range(slover.args['batch_size'])
        )
        
        mask_route = torch.full((x_other.shape[0], 3, 64, 64), -1.0,dtype=torch.float32).to("cpu")
        mask_route[:,1:,:,:] = bg[:,1:,:,:]
        mask_route[:,2:,:,:] = bg[:,2:,:,:]
        for idx, res in enumerate(results):
            mask_route[idx][0],keys = res[0][0],res[1]
            keys_[idx] = keys
        mask_route = mask_route.to(x_other.device)

        # g^* = ▽x_{k+1}( E^o(P(x_{k+1})) - O(P(x_{k+1})))  = (mask_route[:, :, :, :] - x_other[:, :, :, :]) / x_other[:, :, :, :]
        # p_θ(x_k|x_{k+1})p(g*|x_k) = x_other[:, :, :, :] * (mask_route[:, :, :, :] - x_other[:, :, :, :]) / x_other[:, :, :, :] = mask_route[:, :, :, :] - x_other[:, :, :, :]
        # Due to boundary condition, here exp is not used. Meantime, to make the generation process smoother, the (mask_route - x_other) is divided by num_timesteps then added to x_other
        grad = (mask_route[:, :, :, :] - x_other[:, :, :, :]) * t[:, None, None, None] / self.num_timesteps
        x_other[:, :, :, :] = x_other[:, :, :, :] + grad
        if self.args.reparam == 'term':
            if self.direction == 'b':
                coef_t, coef_t_next = self.get_coef_ts(x, t, 1)
                x_1 = x_other
                x_0 = (x - coef_t['coef1'] * x_1) / coef_t['coef0']
            elif self.direction == 'f':
                coef_t, coef_t_next = self.get_coef_ts(x, self.num_timesteps - t, -1)
                x_0 = x_other
                x_1 = (x - coef_t['coef0'] * x_0) / coef_t['coef1']
            x = coef_t_next['coef0'] * x_0 + coef_t_next['coef1'] * x_1
        else:
            x = x_other
        return x,keys_

    def guidance_foward_abl(self,x,t,c_raw,bg,count_,slover: Slover):
        # keys_ = {}
        x_other = self.forward(x, t)
        points_ = (c_raw >= 5 / 6).nonzero(as_tuple=True)
        # print(f'c_raw:{c_raw.shape}')
        x_other[:, 0, :, :][points_] = 1
        if self.args.reparam == 'term':
            if self.direction == 'b':
                coef_t, coef_t_next = self.get_coef_ts(x, t, 1)
                # backward模型的输出,预测x1
                x_1 = x_other
                # x = x0
                # x_0随着时间增大而减小
                x_0 = (x - coef_t['coef1'] * x_1) / coef_t['coef0']
            elif self.direction == 'f':
                coef_t, coef_t_next = self.get_coef_ts(x, self.num_timesteps - t, -1)
                # forward模型的输出，预测x0
                x_0 = x_other
                # x = x1
                x_1 = (x - coef_t['coef0'] * x_0) / coef_t['coef1']
            x = coef_t_next['coef0'] * x_0 + coef_t_next['coef1'] * x_1
        else:
            x = x_other
        return x

    def guidance_inference_abl(self,x,sample=False):
        ones = torch.ones(size=(x.shape[0],), dtype=torch.int64, device=self.device)
        x_raw = x.clone()
        bg = x_raw.clone()
        c_raw = x_raw[:, 0, :, :]
        # 512*1*64*64
        c_raw = (c_raw >= 5 / 6).to(torch.float32)
        x_cache, gt_cache, t_cache= [], [], []
        slover = Slover(args={'batch_size': x.shape[0],'mode':'of'})

        with torch.no_grad():
            count_ = 0
            for t in range(self.num_timesteps):
                tt = ones * t
                x_old = x.clone()
                with torch.autocast(device_type="cuda", enabled=self.args.use_amp):
                    t_old = self.guidance_foward_abl(x, tt, c_raw,bg,count_, slover)
                    # print(f'count_:{count_} done')
                t_old = t_old.float()
                if sample and t == self.num_timesteps - 1:
                    x = t_old
                else:
                    x = t_old + torch.sqrt(2 * self.gammas[t]) * torch.randn_like(x)
                x_cache.append(x.clone())
                if self.args.simplify:
                    if self.args.reparam == 'flow':
                        gt_cache.append((x_raw - x) / (t + 1) * self.num_timesteps)
                    elif self.args.reparam == 'term':
                        gt_cache.append(x_raw)
                    else:
                        gt_cache.append(x_old)
                else:
                    t_new = self._forward(x, tt)
                    gt_cache.append(x + t_old - t_new)
                t_cache.append(self.num_timesteps - 1 - tt)
                count_+=1
        x_cache = torch.stack([x_raw] + x_cache, dim=0).cpu() if sample else torch.cat(x_cache, dim=0).cpu()
        gt_cache = torch.cat(gt_cache, dim=0).cpu()
        t_cache = torch.cat(t_cache, dim=0).cpu()
        return x_cache, gt_cache, t_cache

    def guidance_inferencev1(self,x,mode,alpha,beta,sample=False):
        ones = torch.ones(size=(x.shape[0],), dtype=torch.int64, device=self.device)
        x_raw = x.clone()
        bg = x_raw.clone()
        c_raw = x_raw[:, 0, :, :]
        # 512*1*64*64
        c_raw = (c_raw >= 5 / 6).to(torch.float32)
        x_cache, gt_cache, t_cache, edges_cache = [], [], [], []
        slover = Slover(args={'batch_size': x.shape[0],'mode':mode,'alpha':alpha,'beta':beta})

        with torch.no_grad():
            count_ = 0
            for t in range(self.num_timesteps):
                tt = ones * t
                x_old = x.clone()
                with torch.autocast(device_type="cuda", enabled=self.args.use_amp):
                    t_old,keys = self.guidance_fowardv1(x, tt, c_raw,bg,count_, slover)
                t_old = t_old.float()
                if sample and t == self.num_timesteps - 1:
                    x = t_old
                else:
                    x = t_old + torch.sqrt(2 * self.gammas[t]) * torch.randn_like(x)
                x_cache.append(x.clone())
                edges_cache.append(keys)
                if self.args.simplify:
                    if self.args.reparam == 'flow':
                        gt_cache.append((x_raw - x) / (t + 1) * self.num_timesteps)
                    elif self.args.reparam == 'term':
                        gt_cache.append(x_raw)
                    else:
                        gt_cache.append(x_old)
                else:
                    t_new = self._forward(x, tt)
                    gt_cache.append(x + t_old - t_new)
                t_cache.append(self.num_timesteps - 1 - tt)
                count_+=1
        x_cache = torch.stack([x_raw] + x_cache, dim=0).cpu() if sample else torch.cat(x_cache, dim=0).cpu()
        gt_cache = torch.cat(gt_cache, dim=0).cpu()
        t_cache = torch.cat(t_cache, dim=0).cpu()
        return x_cache, gt_cache, t_cache, edges_cache

    def inference(self, x, sample=False):
        ones = torch.ones(size=(x.shape[0],), dtype=torch.int64, device=self.device)
        x_cache, gt_cache, t_cache = [], [], []
        x_raw = x.clone()
        with torch.no_grad():
            for t in range(self.num_timesteps):
                tt = ones * t
                x_old = x.clone()
                with torch.autocast(device_type="cuda", enabled=self.args.use_amp):
                    t_old = self._forward(x, tt)
                t_old = t_old.float()
                if sample and t == self.num_timesteps - 1:
                    x = t_old
                else:
                    x = t_old + torch.sqrt(2 * self.gammas[t]) * torch.randn_like(x)
                x_cache.append(x.clone())
                if self.args.simplify:
                    if self.args.reparam == 'flow':
                        gt_cache.append((x_raw - x) / (t + 1) * self.num_timesteps)
                    elif self.args.reparam == 'term':
                        gt_cache.append(x_raw)
                    else:
                        gt_cache.append(x_old)
                else:
                    t_new = self._forward(x, tt)
                    gt_cache.append(x + t_old - t_new)
                t_cache.append(self.num_timesteps - 1 - tt)
        x_cache = torch.stack([x_raw] + x_cache, dim=0).cpu() if sample else torch.cat(x_cache, dim=0).cpu()
        gt_cache = torch.cat(gt_cache, dim=0).cpu()
        t_cache = torch.cat(t_cache, dim=0).cpu()
        return x_cache, gt_cache, t_cache


def create_model(name, *args, **kwargs):
    name = name.lower()
    if 'diffusion' in name:
        model = Diffusion
    elif 'flow' in name:
        model = Flow
    elif 'dsb' in name:
        model = DSB
    else:
        raise NotImplementedError
    model = model(*args, **kwargs)
    return model
