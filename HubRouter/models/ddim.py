import torch
import torch.nn as nn
import random
import numpy as np
import lightning.pytorch as pl
# import functools

from utils import ddpm_schedules, instantiate_from_config
# from models import ResnetBlock

# class ContextResnet(pl.LightningModule):
#     def __init__(self, 
#                  input_nc, 
#                  output_nc, 
#                  condition_nc, 
#                  ngf=64, 
#                  n_downsampling=2, 
#                  norm_layer=nn.BatchNorm2d, 
#                  use_dropout=True,
#                  n_blocks=6,
#                  padding_type='zero', 
#                  size=64
#                  ):
#         super().__init__()
#         if type(norm_layer) == functools.partial:
#             use_bias = norm_layer.func == nn.InstanceNorm2d
#         else:
#             use_bias = norm_layer == nn.InstanceNorm2d

#         self.mask_value = torch.zeros([condition_nc, size, size])
#         model = [nn.ZeroPad2d(1),
#                  nn.Conv2d(input_nc+condition_nc, ngf, kernel_size=3, padding=0, bias=use_bias),
#                  norm_layer(ngf),
#                  nn.ReLU(True)]

#         # n_downsampling = 2
#         for i in range(n_downsampling):  # add downsampling layers
#             mult = 2 ** i
#             model += [nn.Conv2d(ngf * mult, ngf * mult * 2, kernel_size=3, stride=2, padding=1, bias=use_bias),
#                       norm_layer(ngf * mult * 2),
#                       nn.ReLU(True)]

#         mult = 2 ** n_downsampling
#         for i in range(n_blocks):  # add ResNet blocks

#             model += [ResnetBlock(ngf * mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
#                                   use_bias=use_bias)]

#         for i in range(n_downsampling):  # add upsampling layers
#             mult = 2 ** (n_downsampling - i)
#             model += [nn.ConvTranspose2d(ngf * mult, int(ngf * mult / 2),
#                                          kernel_size=3, stride=2,
#                                          padding=1, output_padding=1,
#                                          bias=use_bias),
#                       norm_layer(int(ngf * mult / 2)),
#                       nn.ReLU(True)]

#         model += [nn.ZeroPad2d(1)]
#         model += [nn.Conv2d(ngf, output_nc, kernel_size=3, padding=0)]

#         model += [nn.Tanh()]

#         self.model = nn.Sequential(*model)

#     def forward(self, x, c, t, context_mask=None):
#         """Standard forward"""
#         if context_mask is not None:
#             for index in context_mask:
#                 c[index] *= self.mask_value.to(c.device)
#         t = t.view(x.shape[0], 1, 1, 1)
#         x = torch.cat((x, c), 1) + t
#         out = self.model(x)

#         return out

class Diffusion(pl.LightningModule):
    def __init__(self, 
                 unet_config, 
                 betas, 
                 n_T, 
                 ddim_timesteps, 
                 ddim_discr_method, 
                 drop_prob=0.1):
        super().__init__()
        self.save_hyperparameters()
        betas = eval(betas)
        self.nn_model = instantiate_from_config(unet_config)
        self.size = unet_config.params.size
        self.in_channels = unet_config.params.in_channels
        self.ddim_timesteps = ddim_timesteps
        self.ddim_discr_method = ddim_discr_method
        
        self.loss_mse = nn.MSELoss()
        # register_buffer allows accessing dictionary produced by ddpm_schedules
        # e.g. can access self.sqrtab later
        for k, v in ddpm_schedules(betas[0], betas[1], n_T).items():
            self.register_buffer(k, v)

        self.n_T = n_T
        self.drop_prob = drop_prob

    def forward(self, x, c):
        """
        this method is used in training, so samples t and noise randomly
        """

        _ts = torch.randint(1, self.n_T+1, (x.shape[0],)).to(self.device)  # t ~ Uniform(0, n_T)
        
        noise = torch.randn_like(x).to(self.device)  # eps ~ N(0, 1)
        
        x_t = (
            self.sqrtab[_ts, None, None, None] * x
            + self.sqrtmab[_ts, None, None, None] * noise
        )  # This is the x_t, which is sqrt(alphabar) x_0 + sqrt(1-alphabar) * eps
        # We should predict the "error term" from this x_t. Loss is what we return.

        # dropout capacity with some probability
        context_mask = random.choices(range(c.shape[0]), k=int(c.shape[0] * self.drop_prob))
        
        # return MSE between added noise, and our predicted noise
        out = self.nn_model(x_t, c, _ts / self.n_T, context_mask)
        return self.loss_mse(noise, out)
        
    def _extract(self, arr, timesteps, broadcast_shape):
        """
        Extract values from a 1-D numpy array for a batch of indices.
        :param arr: the 1-D numpy array.
        :param timesteps: a tensor of indices into the array to extract.
        :param broadcast_shape: a larger shape of K dimensions with the batch
                                dimension equal to the length of timesteps.
        :return: a tensor of shape [batch_size, 1, ...] where the shape has K dims.
        """
        res = arr.to(device=self.device)[timesteps].float()
        while len(res.shape) < len(broadcast_shape):
            res = res[..., None]
        return res.expand(broadcast_shape)
    
    def ddim_sample(
        self,
        c_sample, 
        ddim_eta=0.0,
        clip_denoised=False, 
        guide_w=0.0):
        
        n_sample = c_sample.shape[0]
        c_i = c_sample.repeat(2,1,1,1)
        # make ddim timestep sequence
        if self.ddim_discr_method == 'uniform':
            interval = self.n_T // self.ddim_timesteps
            ddim_timestep_seq = np.asarray(list(range(1, self.n_T+1, interval)))
        elif self.ddim_discr_method == 'quad':
            ddim_timestep_seq = (
                (np.linspace(0, np.sqrt(self.n_T * .8), self.ddim_timesteps)) ** 2
            ).astype(int)
        else:
            raise NotImplementedError(f'There is no ddim discretization method called "{self.ddim_discr_methodddim_discr_method}"')
        # add one to get the final alpha values right (the ones from first scale to data during sampling)
        ddim_timestep_seq = ddim_timestep_seq + 1
        # previous sequence
        ddim_timestep_prev_seq = np.append(np.array([1]), ddim_timestep_seq[:-1])
        
        # start from pure noise (for each example in the batch)
        sample_img = torch.randn(n_sample, *(self.in_channels, self.size, self.size), device=self.device)

        for i in reversed(range(0, self.ddim_timesteps)):
            # print(f'sampling timestep {i}',end='\r')
            t = torch.full((n_sample,), ddim_timestep_seq[i], device=self.device, dtype=torch.long)
            prev_t = torch.full((n_sample,), ddim_timestep_prev_seq[i], device=self.device, dtype=torch.long)
            
            # 1. get current and previous alpha_cumprod
            alpha_cumprod_t = self._extract(self.alphabar_t, t, sample_img.shape)
            alpha_cumprod_t_prev = self._extract(self.alphabar_t, prev_t, sample_img.shape)
            
            t = t / self.n_T
            prev_t = prev_t / self.n_T
            
            # 2. predict noise using model
            x_i = sample_img.repeat(2,1,1,1)
            t = t.repeat(2,1,1,1)
            eps = self.nn_model(x_i, c_i, t, None)
            eps1 = eps[:n_sample]
            eps2 = eps[n_sample:]
            eps = (1+guide_w)*eps1 - guide_w*eps2
            sample_img = x_i[:n_sample]
            
            # 3. get the predicted x_0
            pred_x0 = (sample_img - torch.sqrt((1. - alpha_cumprod_t)) * eps) / torch.sqrt(alpha_cumprod_t)
            
            if clip_denoised:
                pred_x0 = torch.clamp(pred_x0, min=-1., max=1.)
            
            # 4. compute variance: "sigma_t(η)" -> see formula (16)
            # σ_t = sqrt((1 − α_t−1)/(1 − α_t)) * sqrt(1 − α_t/α_t−1)
            sigmas_t = ddim_eta * torch.sqrt(
                (1 - alpha_cumprod_t_prev) / (1 - alpha_cumprod_t) * (1 - alpha_cumprod_t / alpha_cumprod_t_prev))
            # 5. compute "direction pointing to x_t" of formula (12)
            pred_dir_xt = torch.sqrt(1 - alpha_cumprod_t_prev - sigmas_t**2) * eps
            
            # 6. compute x_{t-1} of formula (12)
            x_prev = torch.sqrt(alpha_cumprod_t_prev) * pred_x0 + pred_dir_xt + sigmas_t * torch.randn_like(sample_img)
            sample_img = x_prev

        return sample_img
    