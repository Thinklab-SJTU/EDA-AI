import os
import numpy as np
import torch
import torch.nn as nn
import torch.optim.lr_scheduler as lrs
import lightning.pytorch as pl
from torchvision.utils import save_image, make_grid
from REST_tool.REST_utils import Evaluator ### REST

from utils import instantiate_from_config, FocalLoss, GANLoss

class DDPMInterface(pl.LightningModule):
    def __init__(self, 
                model_config, 
                model_name, 
                result_dir, 
                ckpt_dir, 
                batch_size, 
                seed):
        super().__init__()
        self.save_hyperparameters()
        self.seed = str(seed)
        
        self.model = instantiate_from_config(model_config)
        self.model.to(self.device)
        
        self.batch_size = batch_size
        self.model_name = model_name
        self.lr = model_config.lr
        self.lr_decay_steps = model_config.lr_decay_steps
        self.lr_decay_rate = model_config.lr_decay_rate
        self.loss_ema = model_config.loss_ema
        self.loss_m = None
        self.samples = True
        self.num_sample = model_config.num_sample
        self.result_dir = result_dir
        self.ckpt_dir = ckpt_dir
        self.ws_test = [float(x) for x in model_config.ws_test.split(',')]
        
    def forward(self, x, c):
        return self.model(x, c)

    def training_step(self, batch, batch_idx):
        self.model.train()
        c = batch['condition'].to(self.device) ### placement and pins
        x = batch['route'].to(self.device) ### routing
    
        loss = self(x, c)
        if self.samples:
            self.c_sample = torch.cat([c[:self.num_sample, :, :, :]])
            self.x_sample = torch.cat([x[:self.num_sample, :, :, :]])
            self.samples = False
            
        self.log_util(loss, 'train_loss')
        return loss
        
    def validation_step(self, batch, batch_idx):
        c = batch['condition'].to(self.device) ### placement and pins
        x = batch['route'].to(self.device) ### routing

        loss = self(x, c)
            
        if self.samples:
            self.c_sample = torch.cat([c[:self.num_sample, :, :, :]])
            self.x_sample = torch.cat([x[:self.num_sample, :, :, :]])
            self.samples = False
            
        self.log_util(loss, 'val_loss')
        return loss

    def test_step(self, batch, batch_idx):
        # Here we just reuse the validation_step for testing
        return self.validation_step(batch, batch_idx)

    def log_util(self, loss, name='loss'):
        self.log(name, loss, on_step=True, on_epoch=True, prog_bar=True, 
            batch_size=self.batch_size)

    def on_validation_epoch_end(self):
        self.samples = True
        self.model.eval()
        print()
        with torch.no_grad():
            for w_i, w in enumerate(self.ws_test):
                x_gen = self.model.ddim_sample(self.c_sample, guide_w=w)
                
                x_all = torch.cat([x_gen, self.x_sample])
                grid = make_grid(x_all*-1 + 1, nrow=self.num_sample)
                
                result_path = os.path.join(self.result_dir, self.model_name + '_seed' + self.seed)
                if not os.path.exists(result_path):
                    os.mkdir(result_path)
                save_image(grid, os.path.join(result_path, f"image_ep{self.current_epoch}_w{w}_val.png"))
                print('saved image at ' + os.path.join(result_path, f"image_ep{self.current_epoch}_w{w}_val.png"))
        if self.current_epoch % 5 == 0:
            torch.save(self.model.state_dict(), self.ckpt_dir + f"{self.model_name}_seed{self.seed}/model_{self.current_epoch}.pth")
            print('saved model at ' + self.ckpt_dir + f"{self.model_name}_seed{self.seed}/model_{self.current_epoch}.pth")
        
    def on_train_epoch_end(self):
        self.samples = True
        self.model.eval()
        print()
        with torch.no_grad():
            for w_i, w in enumerate(self.ws_test):
                x_gen = self.model.ddim_sample(self.c_sample, guide_w=w)
                
                x_all = torch.cat([x_gen, self.x_sample])
                grid = make_grid(x_all*-1 + 1, nrow=self.num_sample)
                
                result_path = os.path.join(self.result_dir, self.model_name + '_seed' + self.seed)
                if not os.path.exists(result_path):
                    os.mkdir(result_path)
                save_image(grid, os.path.join(result_path, f"image_ep{self.current_epoch}_w{w}_train.png"))
                print('saved image at ' + os.path.join(result_path, f"image_ep{self.current_epoch}_w{w}_train.png"))
        
    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        scheduler = lrs.StepLR(optimizer,
                               step_size=self.lr_decay_steps,
                               gamma=self.lr_decay_rate)
        
        return [optimizer], [scheduler]


class GANInterface(pl.LightningModule):
    def __init__(self, 
                 model_config, 
                 model_name, 
                 result_dir, 
                 ckpt_dir, 
                 batch_size, 
                 lambda_norm, 
                 lambda_focal, 
                 seed):
        super().__init__()
        self.save_hyperparameters()
        self.seed = str(seed)
        self.automatic_optimization = False
        self.values = dict() # log_dict
        
        self.generator = instantiate_from_config(model_config.generator_config)
        self.generator.to(self.device)
        self.discriminator = instantiate_from_config(model_config.discriminator_config)
        self.discriminator.to(self.device)
        
        self.batch_size = batch_size
        self.lambda_norm = lambda_norm
        self.lambda_focal = lambda_focal
        self.model_name = model_name
        self.lr = model_config.lr
        self.lr_decay_steps = model_config.lr_decay_steps
        self.lr_decay_rate = model_config.lr_decay_rate
        self.loss_ema = model_config.loss_ema
        self.loss_m = None
        self.samples = True
        self.num_sample = model_config.num_sample
        self.result_dir = result_dir
        self.ckpt_dir = ckpt_dir
        self.criterionGAN = GANLoss('vanilla').to(self.device)
        self.criterionMSE = torch.nn.MSELoss()
        self.focalloss = FocalLoss()
        
    def forward(self, c):
        return self.generator(c)

    def training_step(self, batch, batch_idx):
        self.generator.train()
        self.discriminator.train()
        c = batch['condition'].to(self.device) ### placement and pins
        x = batch['route'].to(self.device) ### routing
        
        if self.samples:
            self.c_sample = torch.cat([c[:self.num_sample, :, :, :]])
            self.x_sample = torch.cat([x[:self.num_sample, :, :, :]])
            self.samples = False
        
        opt_generator, opt_discriminator = self.optimizers()

        loss_D = self.d_loss(x, c)
        self.manual_backward(loss_D)  
        opt_discriminator.step()
        opt_discriminator.zero_grad()

        self.log_util(loss_D, 'train_loss_D')
        
        loss_G = self.g_loss(x, c)
        self.manual_backward(loss_G)
        opt_generator.step()
        opt_generator.zero_grad()
        
        self.log_util(loss_G, 'train_loss_G')
            
        return None
    
    def g_loss(self, x, c):
        fake_x = self(c)
        fake_route_with_condition = torch.cat((c, fake_x), 1)
        pred_fake_route = self.discriminator(fake_route_with_condition)
        loss_G_GAN = self.criterionGAN(pred_fake_route, True)
        loss_G_norm = self.criterionMSE(fake_x, x) * self.lambda_norm
        loss_G_focal = self.focalloss((fake_x + 1.) / 2., (x + 1.) / 2.) * self.lambda_focal
        loss_G = (loss_G_norm + loss_G_focal) * 5000 + loss_G_GAN
        return loss_G
        
    def d_loss(self, x, c):
        fake_route_with_condition = torch.cat((c, self(c)), 1)
        pred_fake_route = self.discriminator(fake_route_with_condition.detach())
        loss_D_fake = self.criterionGAN(pred_fake_route, False)
        
        real_route_with_condition = torch.cat((c, x), 1)
        pred_real_route = self.discriminator(real_route_with_condition)
        loss_D_real = self.criterionGAN(pred_real_route, True)
        loss_D = (loss_D_fake + loss_D_real) * 0.5
        return loss_D
    
    def validation_step(self, batch, batch_idx):
        self.generator.eval()
        c = batch['condition'].to(self.device) ### placement and pins
        x = batch['route'].to(self.device) ### routing
            
        if self.samples:
            self.c_sample = torch.cat([c[:self.num_sample, :, :, :]])
            self.x_sample = torch.cat([x[:self.num_sample, :, :, :]])
            self.samples = False

        loss_G = self.g_loss(x, c)
        loss_D = self.d_loss(x, c)
        self.log_util(loss_G, 'val_loss_G')
        self.log_util(loss_D, 'val_loss_D')
        
        return loss_G+loss_D

    def test_step(self, batch, batch_idx):
        # Here we just reuse the validation_step for testing
        return self.validation_step(batch, batch_idx)

    def log_util(self, loss, name='loss'):
        self.values[name] = loss
        self.log_dict(self.values, logger=True, prog_bar=True, on_step=True, on_epoch=True, 
                      batch_size=self.batch_size)

    def on_validation_epoch_end(self):
        self.samples = True
        self.generator.eval()
        print()
        with torch.no_grad():
            x_gen = self(self.c_sample)
            x_all = torch.cat([x_gen, self.x_sample])
            grid = make_grid(x_all*-1 + 1, nrow=self.num_sample)
                
            result_path = os.path.join(self.result_dir, self.model_name + '_seed' + self.seed)
            if not os.path.exists(result_path):
                os.mkdir(result_path)
            save_image(grid, os.path.join(result_path, f"image_ep{self.current_epoch}_val.png"))
            print('saved image at ' + os.path.join(result_path, f"image_ep{self.current_epoch}_val.png"))
            
        if self.current_epoch % 5 == 0:
            self.discriminator.eval()
            torch.save(self.generator.state_dict(), self.ckpt_dir + f"{self.model_name}_seed{self.seed}/generator_{self.current_epoch}_seed{self.seed}.pth")
            torch.save(self.discriminator.state_dict(), self.ckpt_dir + f"{self.model_name}_seed{self.seed}/discriminator_{self.current_epoch}.pth")
            print('saved model at ' + self.ckpt_dir + f"{self.model_name}_seed{self.seed}/generator_{self.current_epoch}.pth")
            print('saved model at ' + self.ckpt_dir + f"{self.model_name}_seed{self.seed}/discriminator_{self.current_epoch}.pth")
        
    def on_train_epoch_end(self):
        self.samples = True
        self.generator.eval()
        print()
        with torch.no_grad():
            x_gen = self(self.c_sample)
            x_all = torch.cat([x_gen, self.x_sample])
            grid = make_grid(x_all*-1 + 1, nrow=self.num_sample)
            
            result_path = os.path.join(self.result_dir, self.model_name + '_seed' + self.seed)
            if not os.path.exists(result_path):
                os.mkdir(result_path)
            save_image(grid, os.path.join(result_path, f"image_ep{self.current_epoch}_train.png"))
            print('saved image at ' + os.path.join(result_path, f"image_ep{self.current_epoch}_train.png"))
        
    def configure_optimizers(self):
        opt_generator = torch.optim.Adam(self.generator.parameters(), lr=self.lr, betas=(0.5, 0.999),
                                                weight_decay=0.01)
        opt_discriminator = torch.optim.Adam(self.discriminator.parameters(), lr=self.lr, betas=(0.5, 0.999), weight_decay=0.01)
        scheduler_generator = lrs.StepLR(opt_generator,
                                         step_size=self.lr_decay_steps,
                                         gamma=self.lr_decay_rate)
        scheduler_discriminator = lrs.StepLR(opt_discriminator,
                                             step_size=self.lr_decay_steps,
                                             gamma=self.lr_decay_rate)
                                         
        return [opt_generator, opt_discriminator], [scheduler_generator, scheduler_discriminator]
        
class VAEInterface(pl.LightningModule):
    def __init__(self, 
                 model_config, 
                 model_name, 
                 result_dir, 
                 ckpt_dir, 
                 batch_size, 
                 seed, 
                 latent_size):
        super().__init__()
        self.save_hyperparameters()
        self.seed = str(seed)
        self.values = dict() # log_dict
        
        self.vae = instantiate_from_config(model_config.vae_config)
        self.vae.to(self.device)
        self.focalloss = FocalLoss()
        
        self.batch_size = batch_size
        self.latent_size = latent_size
        self.model_name = model_name
        self.lr = model_config.lr
        self.lr_decay_steps = model_config.lr_decay_steps
        self.lr_decay_rate = model_config.lr_decay_rate
        self.loss_ema = model_config.loss_ema
        self.loss_m = None
        self.samples = True
        self.num_sample = model_config.num_sample
        self.result_dir = result_dir
        self.ckpt_dir = ckpt_dir
        self.criterionMSE = torch.nn.MSELoss()
            
    def forward(self, x, c):
        means, log_var = self.vae(x, c)
        z = self.reparameterize(means, log_var)
        recon_x = self.vae.inference(z, c)
        return recon_x, means, log_var, z
        
    def loss_fn(self, x, c):
        recon_x, mean, log_var, z = self(x, c)
        BCE = self.criterionMSE(recon_x.view(x.size(0), -1), x.view(x.size(0), -1)) * 100
        focal = self.focalloss((recon_x.view(x.size(0), -1) + 1.) / 2., (x.view(x.size(0), -1) + 1.) / 2.) * 100
        KLD = -0.5 * torch.sum(1 + log_var - mean.pow(2) - log_var.exp())

        return ((BCE + focal)*5000 + KLD) / x.size(0)
        
    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + eps * std

    def training_step(self, batch, batch_idx):
        self.vae.train()
        c = batch['condition'].to(self.device) ### placement and pins
        x = batch['route'].to(self.device) ### routing
        
        if self.samples:
            self.c_sample = torch.cat([c[:self.num_sample, :, :, :]])
            self.x_sample = torch.cat([x[:self.num_sample, :, :, :]])
            self.samples = False
        
        loss = self.loss_fn(x, c)
        self.log_util(loss, 'train_loss')
        return loss
    
    def validation_step(self, batch, batch_idx):
        self.vae.eval()
        c = batch['condition'].to(self.device) ### placement and pins
        x = batch['route'].to(self.device) ### routing
            
        if self.samples:
            self.c_sample = torch.cat([c[:self.num_sample, :, :, :]])
            self.x_sample = torch.cat([x[:self.num_sample, :, :, :]])
            self.samples = False

        loss = self.loss_fn(x, c)
        self.log_util(loss, 'val_loss')
        return loss

    def test_step(self, batch, batch_idx):
        # Here we just reuse the validation_step for testing
        return self.validation_step(batch, batch_idx)

    def log_util(self, loss, name='loss'):
        self.values[name] = loss
        self.log_dict(self.values, logger=True, prog_bar=True, on_step=True, on_epoch=True, 
                      batch_size=self.batch_size)

    def on_validation_epoch_end(self):
        self.samples = True
        self.vae.eval()
        print()
        with torch.no_grad():
            z = torch.randn([self.c_sample.size(0), self.latent_size]).to(self.device)
            x_gen = self.vae.inference(z, c=self.c_sample)
            x_all = torch.cat([x_gen, self.x_sample])
            grid = make_grid(x_all*-1 + 1, nrow=self.num_sample)
                
            result_path = os.path.join(self.result_dir, self.model_name + '_seed' + self.seed)
            if not os.path.exists(result_path):
                os.mkdir(result_path)
            save_image(grid, os.path.join(result_path, f"image_ep{self.current_epoch}_val.png"))
            print('saved image at ' + os.path.join(result_path, f"image_ep{self.current_epoch}_val.png"))
            
        if self.current_epoch % 5 == 0:
            torch.save(self.vae.state_dict(), self.ckpt_dir + f"{self.model_name}_seed{self.seed}/vae_{self.current_epoch}_seed{self.seed}.pth")
            print('saved model at ' + self.ckpt_dir + f"{self.model_name}_seed{self.seed}/vae_{self.current_epoch}.pth")
        
    def on_train_epoch_end(self):
        self.samples = True
        self.vae.eval()
        print()
        with torch.no_grad():
            z = torch.randn([self.c_sample.size(0), self.latent_size]).to(self.device)
            x_gen = self.vae.inference(z, c=self.c_sample)
            x_all = torch.cat([x_gen, self.x_sample])
            grid = make_grid(x_all*-1 + 1, nrow=self.num_sample)
            
            result_path = os.path.join(self.result_dir, self.model_name + '_seed' + self.seed)
            if not os.path.exists(result_path):
                os.mkdir(result_path)
            save_image(grid, os.path.join(result_path, f"image_ep{self.current_epoch}_train.png"))
            print('saved image at ' + os.path.join(result_path, f"image_ep{self.current_epoch}_train.png"))
        
    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=self.lr)
        scheduler = lrs.StepLR(optimizer,
                               step_size=self.lr_decay_steps,
                               gamma=self.lr_decay_rate)
        
        return [optimizer], [scheduler]
        
class ACInterface(pl.LightningModule):
    def __init__(self, 
                 model_config, 
                 model_name, 
                 result_dir, 
                 ckpt_dir, 
                 batch_size, 
                 degree, 
                 seed):
        super().__init__()
        self.save_hyperparameters()
        self.seed = str(seed)
        self.values = dict() # log_dict
        
        self.actor = instantiate_from_config(model_config.actor_config)
        self.actor.to(self.device)
        self.critic = instantiate_from_config(model_config.critic_config)
        self.critic.to(self.device)
        
        
        self.batch_size = batch_size
        self.model_name = model_name
        self.lr = model_config.lr
        self.degree = degree
        self.result_dir = result_dir
        self.ckpt_dir = ckpt_dir
        self.evaluator = Evaluator(path=os.getcwd() + '/REST_tool/algorithms/libeval.so') # We borrow the evaluator from REST
        self.mse_loss = torch.nn.MSELoss()
        
    def forward(self, input, trans=False):
        outputs, log_probs = self.actor(input, trans)
        return outputs, log_probs

    def training_step(self, batch, batch_idx):
        self.actor.train()
        self.critic.train()
        input = batch.to(self.device)
        outputs, log_probs = self(input)
        predictions = self.critic(input)
    
        lengths = self.evaluator.eval_batch(batch.cpu().detach().numpy(), outputs.cpu().detach().numpy(), self.degree)
        length_tensor = lengths.detach().type(torch.float).to(self.device)
        
        with torch.no_grad():
            disadvantage = length_tensor - predictions
        actor_loss = torch.mean(disadvantage * log_probs)
        critic_loss = self.mse_loss(predictions, length_tensor)
        loss = actor_loss + critic_loss
        self.log_util(loss, 'train_loss')
            
        return loss
    
    def validation_step(self, batch, batch_idx):
        self.actor.eval()
        input = batch.to(self.device)
        with torch.no_grad():
            outputs, _ = self(input, True)
        lengths = self.evaluator.eval_batch(batch.cpu().detach().numpy(), outputs.cpu().detach().numpy(), self.degree)
        self.log_util(lengths.mean(), 'val_lengths')
        
        return lengths

    def test_step(self, batch, batch_idx):
        # Here we just reuse the validation_step for testing
        return self.validation_step(batch, batch_idx)

    def log_util(self, loss, name='loss'):
        self.values[name] = loss
        self.log_dict(self.values, logger=True, prog_bar=True, on_step=True, on_epoch=True, 
                      batch_size=self.batch_size)
        
    def on_train_epoch_end(self):
        torch.nn.utils.clip_grad_norm_(self.actor.parameters(), 1.)
        torch.nn.utils.clip_grad_norm_(self.critic.parameters(), 1.)
        
    def configure_optimizers(self):
        optimizer = torch.optim.Adam(list(self.actor.parameters()) + list(self.critic.parameters()), lr=self.lr, eps=1e-5)
        
        return [optimizer], []