import torch.nn as nn
import lightning.pytorch as pl
import functools

from models import ResnetBlock

class ResnetGenerator(pl.LightningModule):
    """Resnet-based generator that consists of Resnet blocks between a few downsampling/upsampling operations.
    We adapt Torch code and idea from Justin Johnson's neural style transfer project(https://github.com/jcjohnson/fast-neural-style)
    """
    def __init__(self, 
                 input_nc, 
                 output_nc, 
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
            input_nc (int)      -- the number of channels in input images
            output_nc (int)     -- the number of channels in output images
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
        #          nn.Conv2d(input_nc, ngf, kernel_size=7, padding=0, bias=use_bias),
        #          norm_layer(ngf),
        #          nn.ReLU(True)]

        # modified
        model = [nn.ZeroPad2d(1),
                 nn.Conv2d(input_nc, ngf, kernel_size=3, padding=0, bias=use_bias),
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

        for i in range(n_downsampling):  # add upsampling layers
            mult = 2 ** (n_downsampling - i)
            model += [nn.ConvTranspose2d(ngf * mult, int(ngf * mult / 2),
                                         kernel_size=3, stride=2,
                                         padding=1, output_padding=1,
                                         bias=use_bias),
                      norm_layer(int(ngf * mult / 2)),
                      nn.ReLU(True)]
        # original
        # model += [nn.ReflectionPad2d(3)]
        # model += [nn.Conv2d(ngf, output_nc, kernel_size=7, padding=0)]

        # modified
        model += [nn.ZeroPad2d(1)]
        model += [nn.Conv2d(ngf, output_nc, kernel_size=3, padding=0)]

        model += [nn.Tanh()]

        self.model = nn.Sequential(*model)

    def forward(self, input):
        """Standard forward"""
        return self.model(input)
        
class Discriminator(pl.LightningModule):
    def __init__(self, 
                 input_nc, 
                 ndf=64, 
                 n_layers=3, 
                 norm_layer=nn.BatchNorm2d, 
                 use_dropout=False):
        self.save_hyperparameters()
        """Construct a PatchGAN discriminator

        Parameters:
            input_nc (int)  -- the number of channels in input images
            ndf (int)       -- the number of filters in the last conv layer
            n_layers (int)  -- the number of conv layers in the discriminator
            norm_layer      -- normalization layer
        """
        super(Discriminator, self).__init__()
        if type(norm_layer) == functools.partial:  # no need to use bias as BatchNorm2d has affine parameters
            use_bias = norm_layer.func == nn.InstanceNorm2d
        else:
            use_bias = norm_layer == nn.InstanceNorm2d

        kw = 4
        padw = 1
        padding_type = 'zero'

        sequence = [nn.ZeroPad2d(1),
                    nn.Conv2d(input_nc, ndf, kernel_size=3, padding=0, bias=use_bias),
                    norm_layer(ndf),
                    nn.ReLU(True)]
        nf_mult = 1
        nf_mult_prev = 1
        for n in range(1, n_layers):  # gradually increase the number of filters
            nf_mult_prev = nf_mult
            nf_mult = min(2 ** n, 8)
            sequence += [
                nn.Conv2d(ndf * nf_mult_prev, ndf * nf_mult, kernel_size=kw, stride=2, padding=padw, bias=use_bias),
                norm_layer(ndf * nf_mult),
                nn.LeakyReLU(0.2, True),
                ResnetBlock(ndf * nf_mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                            use_bias=use_bias),
                ResnetBlock(ndf * nf_mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                            use_bias=use_bias),
                ResnetBlock(ndf * nf_mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                            use_bias=use_bias)
            ]

        self.preceded_net = nn.Sequential(*sequence)

        nf_mult_prev = nf_mult
        nf_mult = min(2 ** n_layers, 8)

        self.d_net = nn.Sequential(
            nn.Conv2d(ndf * nf_mult_prev, ndf * nf_mult, kernel_size=kw, stride=2, padding=padw, bias=use_bias),
            norm_layer(ndf * nf_mult),
            nn.LeakyReLU(0.2, True),
            ResnetBlock(ndf * nf_mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                        use_bias=use_bias),
            ResnetBlock(ndf * nf_mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                        use_bias=use_bias),
            ResnetBlock(ndf * nf_mult, padding_type=padding_type, norm_layer=norm_layer, use_dropout=use_dropout,
                        use_bias=use_bias),
            nn.Conv2d(ndf * nf_mult, 1, kernel_size=kw, stride=1, padding=1, bias=use_bias)
        )

    def forward(self, input):
        """Standard forward."""
        preceded_fm = self.preceded_net(input)
        out = self.d_net(preceded_fm)
        return out   
        
    