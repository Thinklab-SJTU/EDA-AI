import torch
import torch.nn as nn
import lightning.pytorch as pl
import functools

from models import EmbedFC, ResnetBlock

class VAE(pl.LightningModule):
    def __init__(self, 
                 in_channels, 
                 out_channels, 
                 latent_size=64, 
                 ngf=64, 
                 n_downsampling=2, 
                 norm_layer=nn.BatchNorm2d, 
                 use_dropout=True,
                 n_blocks=6,
                 padding_type='zero', 
                 size=64, 
                 output_nc=1, 
                 c_in_channels=3):
        super().__init__()
        self.save_hyperparameters()
        self.encoder = Encoder(in_channels=in_channels, 
                               latent_size=latent_size, 
                               ngf=ngf, 
                               n_downsampling=n_downsampling, 
                               norm_layer=nn.BatchNorm2d, 
                               use_dropout=use_dropout,
                               n_blocks=n_blocks,
                               padding_type=padding_type).to(self.device)
        self.decoder = Decoder(out_channels=out_channels, 
                               c_in_channels=c_in_channels, 
                               n_downsampling=n_downsampling, 
                               ngf=ngf, 
                               latent_size=latent_size, 
                               size=size).to(self.device)
    
    def forward(self, x, c):
        return self.encoder(x, c)
        
    def inference(self, z, c):
        return self.decoder(z, c)
    
    
class Encoder(pl.LightningModule):
    """Resnet-based generator that consists of Resnet blocks between a few downsampling/upsampling operations.
    We adapt Torch code and idea from Justin Johnson's neural style transfer project(https://github.com/jcjohnson/fast-neural-style)
    """
    def __init__(self, 
                 in_channels, 
                 latent_size=64, 
                 ngf=64, 
                 n_downsampling=2, 
                 norm_layer=nn.BatchNorm2d, 
                 use_dropout=True,
                 n_blocks=6,
                 padding_type='zero'):
        super().__init__()
        self.save_hyperparameters()
        """Construct a Resnet-based generator

        Parameters:
            in_channels (int)      -- the number of channels in input images
            output_nc (int)     -- the number of channels in hidden layers
            ngf (int)           -- the number of filters in the last conv layer
            norm_layer          -- normalization layer
            use_dropout (bool)  -- if use dropout layers
            n_blocks (int)      -- the number of ResNet blocks
            padding_type (str)  -- the name of padding layer in conv layers: reflect | replicate | zero
        """
        if type(norm_layer) == functools.partial:
            use_bias = norm_layer.func == nn.InstanceNorm2d
        else:
            use_bias = norm_layer == nn.InstanceNorm2d

        # original
        # model = [nn.ReflectionPad2d(3),
        #          nn.Conv2d(in_channels, ngf, kernel_size=7, padding=0, bias=use_bias),
        #          norm_layer(ngf),
        #          nn.ReLU(True)]

        # modified
        model = [nn.ZeroPad2d(1),
                 nn.Conv2d(in_channels, ngf, kernel_size=3, padding=0, bias=use_bias),
                 norm_layer(ngf),
                 nn.ReLU(True)]

        # n_downsampling = 2
        for i in range(n_downsampling):  # add downsampling layers
            mult = 2 ** i
            model += [nn.Conv2d(ngf * mult, ngf * mult * 2, kernel_size=3, stride=2, padding=1, bias=use_bias),
                      norm_layer(ngf * mult * 2),
                      nn.ReLU(True)]

        mult = 2 ** n_downsampling
        for i in range(n_blocks):  # add ResNet blocks

            model += [ResnetBlock(ngf * mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                                  use_bias=use_bias)]

        self.model = nn.Sequential(*model)
        
        self.linear_means = EmbedFC(ngf * mult * (ngf//4)**2, 2*ngf, latent_size)
        self.linear_log_var = EmbedFC(ngf * mult * (ngf//4)**2, 2*ngf, latent_size)

    def forward(self, x, c):
        """Standard forward"""
        x = torch.cat((x, c), 1)
        output = self.model(x)
        return self.linear_means(output), self.linear_log_var(output)
        
class Decoder(pl.LightningModule):
    def __init__(self, 
                 out_channels, 
                 c_in_channels, 
                 n_downsampling = 2, 
                 norm_layer=nn.BatchNorm2d, 
                 ngf=64, 
                 latent_size=64, 
                 size=64):
        super().__init__()
        self.save_hyperparameters()
        if type(norm_layer) == functools.partial:
            use_bias = norm_layer.func == nn.InstanceNorm2d
        else:
            use_bias = norm_layer == nn.InstanceNorm2d
        
        self.size = size
        self.ngf = ngf
        self.n_downsampling = n_downsampling
        mult = 2**self.n_downsampling
        self.zemb = EmbedFC(latent_size, 2*ngf, self.ngf * mult//2 * (ngf//4)**2)
        self.cemb = EmbedFC(c_in_channels*size*size, 2*ngf, self.ngf * mult//2 * (ngf//4)**2)
        
        # model = [nn.Conv2d(ngf * mult + c_in_channels, ngf * mult, kernel_size=3, stride=2, padding=1, bias=use_bias),
                 # norm_layer(ngf * mult),
                 # nn.ReLU(True)]
        model = []
        
        for i in range(n_downsampling):  # add upsampling layers
            mult = 2 ** (n_downsampling - i)
            model += [nn.ConvTranspose2d(ngf * mult, int(ngf * mult / 2),
                                         kernel_size=3, stride=2,
                                         padding=1, output_padding=1,
                                         bias=use_bias),
                      norm_layer(int(ngf * mult / 2)),
                      nn.ReLU(True)]
                      
        model += [nn.ZeroPad2d(1)]
        model += [nn.Conv2d(ngf, out_channels, kernel_size=3, padding=0)]
        
        model += [nn.Tanh()]
        
        self.model = nn.Sequential(*model)

    def forward(self, z, c):
        zemb = self.zemb(z).view(z.shape[0], self.ngf * (2**self.n_downsampling) //2, self.ngf//4, self.ngf//4)
        cemb = self.cemb(c).view(c.shape[0], self.ngf * (2**self.n_downsampling) //2, self.ngf//4, self.ngf//4)
        z = torch.cat([zemb,cemb], axis=1)
        out = self.model(z)
        return out