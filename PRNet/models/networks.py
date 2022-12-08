import torch
import torch.nn as nn
from torch.nn import init
import functools
from torch.optim import lr_scheduler


def get_norm_layer(norm_type='instance'):
    """Return a normalization layer

    Parameters:
        norm_type (str) -- the name of the normalization layer: batch | instance | none

    For BatchNorm, we use learnable affine parameters and track running statistics (mean/stddev).
    For InstanceNorm, we do not use learnable affine parameters. We do not track running statistics.
    """
    if norm_type == 'batch':
        norm_layer = functools.partial(nn.BatchNorm2d, affine=True, track_running_stats=True)
    elif norm_type == 'instance':
        norm_layer = functools.partial(nn.InstanceNorm2d, affine=False, track_running_stats=False)
    elif norm_type == 'none':
        def norm_layer(x):
            return Identity()
    else:
        raise NotImplementedError('normalization layer [%s] is not found' % norm_type)
    return norm_layer


def get_scheduler(optimizer, opt):
    """Return a learning rate scheduler

    Parameters:
        optimizer          -- the optimizer of the network
        opt (option class) -- stores all the experiment flags; needs to be a subclass of BaseOptions．　
                              opt.lr_policy is the name of learning rate policy: linear | step | plateau | cosine

    For 'linear', we keep the same learning rate for the first <opt.n_epochs> epochs
    and linearly decay the rate to zero over the next <opt.n_epochs_decay> epochs.
    For other schedulers (step, plateau, and cosine), we use the default PyTorch schedulers.
    See https://pytorch.org/docs/stable/optim.html for more details.
    """
    if opt.lr_policy == 'linear':
        def lambda_rule(epoch):
            lr_l = 1.0 - max(0, epoch + opt.epoch_count - opt.n_epochs) / float(opt.n_epochs_decay + 1)
            return lr_l

        scheduler = lr_scheduler.LambdaLR(optimizer, lr_lambda=lambda_rule)
    elif opt.lr_policy == 'step':
        scheduler = lr_scheduler.StepLR(optimizer, step_size=opt.lr_decay_iters, gamma=0.1)
    elif opt.lr_policy == 'plateau':
        scheduler = lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.2, threshold=0.01, patience=5)
    elif opt.lr_policy == 'cosine':
        scheduler = lr_scheduler.CosineAnnealingLR(optimizer, T_max=opt.n_epochs, eta_min=0)
    else:
        return NotImplementedError('learning rate policy [%s] is not implemented', opt.lr_policy)
    return scheduler


def init_weights(net, init_type='normal', init_gain=0.02):
    """Initialize network weights.

    Parameters:
        net (network)   -- network to be initialized
        init_type (str) -- the name of an initialization method: normal | xavier | kaiming | orthogonal
        init_gain (float)    -- scaling factor for normal, xavier and orthogonal.
    """

    def init_func(m):  # define the initialization function
        classname = m.__class__.__name__
        if hasattr(m, 'weight') and (classname.find('Conv') != -1 or classname.find('Linear') != -1):
            if init_type == 'normal':
                init.normal_(m.weight.data, 0.0, init_gain)
            elif init_type == 'xavier':
                init.xavier_normal_(m.weight.data, gain=init_gain)
            elif init_type == 'kaiming':
                init.kaiming_normal_(m.weight.data, a=0, mode='fan_in')
            elif init_type == 'orthogonal':
                init.orthogonal_(m.weight.data, gain=init_gain)
            else:
                raise NotImplementedError('initialization method [%s] is not implemented' % init_type)
            if hasattr(m, 'bias') and m.bias is not None:
                init.constant_(m.bias.data, 0.0)
        elif classname.find(
                'BatchNorm2d') != -1:  # BatchNorm Layer's weight is not a matrix; only normal distribution applies.
            init.normal_(m.weight.data, 1.0, init_gain)
            init.constant_(m.bias.data, 0.0)

    print('initialize network with %s' % init_type)
    net.apply(init_func)  # apply the initialization function <init_func>


def init_net(net, init_type='normal', init_gain=0.02, gpu_ids=[]):
    """Initialize a network: 1. register CPU/GPU device (with multi-GPU support); 2. initialize the network weights
    Parameters:
        net (network)      -- the network to be initialized
        init_type (str)    -- the name of an initialization method: normal | xavier | kaiming | orthogonal
        gain (float)       -- scaling factor for normal, xavier and orthogonal.
        gpu_ids (int list) -- which GPUs the network runs on: e.g., 0,1,2

    Return an initialized network.
    """
    if len(gpu_ids) > 0:
        assert (torch.cuda.is_available())
        net.to(gpu_ids[0])
        net = torch.nn.DataParallel(net, gpu_ids)  # multi-GPUs
    init_weights(net, init_type, init_gain=init_gain)
    return net


def define_G(input_nc, output_nc, ngf, netG, norm='batch', use_dropout=True, init_type='normal', init_gain=0.02,
             gpu_ids=[]):
    """Create a generator

    Parameters:
        input_nc (int) -- the number of channels in input images
        output_nc (int) -- the number of channels in output images
        ngf (int) -- the number of filters in the last conv layer
        netG (str) -- the architecture's name:   |  |  |
        norm (str) -- the name of normalization layers used in the network: batch | instance | none
        use_dropout (bool) -- if use dropout layers.
        init_type (str)    -- the name of our initialization method.
        init_gain (float)  -- scaling factor for normal, xavier and orthogonal.
        gpu_ids (int list) -- which GPUs the network runs on: e.g., 0,1,2

    Returns a generator

    Our current implementation provides two types of generators:
        Basic (for 64x64 input images) and Expanded (for 128x128 input images)
    The generator has been initialized by <init_net>. It uses RELU for non-linearity.
    """
    net = None
    norm_layer = nn.BatchNorm2d
    if netG == 'g_128':
        net = Filling(input_nc, output_nc)
    elif netG == 'baseline_128':
        net = ResnetGenerator(input_nc, output_nc, ngf=32, n_downsampling=3, norm_layer=norm_layer,
                              use_dropout=use_dropout, n_blocks=12)
    elif netG == 'g_128_f':
        net = SimpleEnhancer(input_nc, output_nc)
    elif netG == 'resnet_9blocks':
        net = ResnetGenerator(input_nc, output_nc, ngf, norm_layer=norm_layer, use_dropout=use_dropout, n_blocks=9)
    elif netG == 'cvae':
        net = CVAE(input_nc, output_nc, ngf, norm_layer=norm_layer, use_dropout=use_dropout)
    elif netG == 'resnet_6blocks':
        net = ResnetGenerator(input_nc, output_nc, ngf, norm_layer=norm_layer, use_dropout=use_dropout, n_blocks=6)
    elif netG == 'unet_64':
        net = UnetGenerator(input_nc, output_nc, 6, ngf, norm_layer=norm_layer,
                            use_dropout=use_dropout)
    # elif netG == 'unet_64_sa':
    #     net = UnetGenerator(input_nc, output_nc, 6, ngf, norm_layer=norm_layer,
    #                         use_dropout=use_dropout, attention=True)
    else:
        raise NotImplementedError('Generator model name [%s] is not recognized' % netG)
    return init_net(net, init_type, init_gain, gpu_ids)


def define_D(input_nc, ndf, netD, n_layers_D=5, norm='batch', init_type='normal', init_gain=0.02, gpu_ids=[]):
    """Create a discriminator.

    Parameters:
        input_nc (int)     -- the number of channels in input images
        ndf (int)          -- the number of filters in the first conv layer
        netD (str)         -- the architecture's name: basic | n_layers | pixel
        n_layers_D (int)   -- the number of conv layers in the discriminator; effective when netD=='n_layers'
        norm (str)         -- the type of normalization layers used in the network.
        init_type (str)    -- the name of the initialization method.
        init_gain (float)  -- scaling factor for normal, xavier and orthogonal.
        gpu_ids (int list) -- which GPUs the network runs on: e.g., 0,1,2

    Returns a discriminator

    Our current implementation provides three types of discriminators:
        [basic]: 'PatchGAN' classifier described in the original pix2pix paper.
        Such a patch-level discriminator architecture has fewer parameters
        than a full-image discriminator and can work on arbitrarily-sized images
        in a fully convolutional fashion.

        [critic]: This classifier judges the output images merely by the connectivity of the generated image.

        [multi]: This classifier extends the basic one with a critic classifier that discriminates the output images
        from both realness and connectivity.

    The discriminator has been initialized by <init_net>. It uses Leakly RELU for non-linearity.
    """
    net = None
    norm_layer = nn.BatchNorm2d

    if netD == 'basic':  # default classifier
        net = Discriminator(input_nc, ndf, n_layers=3, norm_layer=norm_layer)
    elif netD == 'critic':
        net = Critic(input_nc, ndf, n_layers=3, norm_layer=norm_layer)
    elif netD == 'multi':
        net = Multi_Discriminator(input_nc, ndf, n_layers=3, norm_layer=norm_layer)
    else:
        raise NotImplementedError('Discriminator model name [%s] is not recognized' % netD)
    return init_net(net, init_type, init_gain, gpu_ids)


##############################################################################
# Classes
##############################################################################
class GANLoss(nn.Module):
    """Define different GAN objectives.

    The GANLoss class abstracts away the need to create the target label tensor
    that has the same size as the input.
    """

    def __init__(self, gan_mode, target_real_label=1., target_fake_label=0.):
        """ Initialize the GANLoss class.

        Parameters:
            gan_mode (str) - - the type of GAN objective. It currently supports vanilla, lsgan, and wgangp.
            target_real_label (bool) - - label for a real image
            target_fake_label (bool) - - label of a fake image

        Note: Do not use sigmoid as the last layer of Discriminator.
        LSGAN needs no sigmoid. vanilla GANs will handle it with BCEWithLogitsLoss.
        """
        super(GANLoss, self).__init__()
        self.register_buffer('real_label', torch.tensor(target_real_label))
        self.register_buffer('fake_label', torch.tensor(target_fake_label))
        self.gan_mode = gan_mode
        if gan_mode == 'lsgan':
            self.loss = nn.L1Loss()
        elif gan_mode == 'vanilla':
            self.loss = nn.BCEWithLogitsLoss()
        elif gan_mode in ['wgangp']:
            self.loss = None
        else:
            raise NotImplementedError('gan mode %s not implemented' % gan_mode)

    def get_target_tensor(self, prediction, target_is_real):
        """Create label tensors with the same size as the input.

        Parameters:
            prediction (tensor) - - typically the prediction from a discriminator
            target_is_real (bool) - - if the ground truth label is for real images or fake images

        Returns:
            A label tensor filled with ground truth label, and with the size of the input
        """

        if target_is_real:
            target_tensor = self.real_label
        else:
            target_tensor = self.fake_label
        return target_tensor.expand_as(prediction)

    def __call__(self, prediction, target_is_real):
        """Calculate loss given Discriminator's output and ground truth labels.

        Parameters:
            prediction (tensor) - - typically the prediction output from a discriminator
            target_is_real (bool) - - if the ground truth label is for real images or fake images

        Returns:
            the calculated loss.
        """
        if self.gan_mode in ['lsgan', 'vanilla']:
            target_tensor = self.get_target_tensor(prediction, target_is_real)
            loss = self.loss(prediction, target_tensor)
        elif self.gan_mode == 'wgangp':
            if target_is_real:
                loss = -prediction.mean()
            else:
                loss = prediction.mean()
        return loss


def cal_gradient_penalty(netD, real_data, fake_data, device, type='mixed', constant=1.0, lambda_gp=10.0):
    """Calculate the gradient penalty loss, used in WGAN-GP paper https://arxiv.org/abs/1704.00028

    Arguments:
        netD (network)              -- discriminator network
        real_data (tensor array)    -- real images
        fake_data (tensor array)    -- generated images from the generator
        device (str)                -- GPU / CPU: from torch.device('cuda:{}'.format(self.gpu_ids[0])) if self.gpu_ids else torch.device('cpu')
        type (str)                  -- if we mix real and fake data or not [real | fake | mixed].
        constant (float)            -- the constant used in formula ( ||gradient||_2 - constant)^2
        lambda_gp (float)           -- weight for this loss

    Returns the gradient penalty loss
    """
    if lambda_gp > 0.0:
        if type == 'real':  # either use real images, fake images, or a linear interpolation of two.
            interpolatesv = real_data
        elif type == 'fake':
            interpolatesv = fake_data
        elif type == 'mixed':
            alpha = torch.rand(real_data.shape[0], 1, device=device)
            alpha = alpha.expand(real_data.shape[0], real_data.nelement() // real_data.shape[0]).contiguous().view(
                *real_data.shape)
            interpolatesv = alpha * real_data + ((1 - alpha) * fake_data)
        else:
            raise NotImplementedError('{} not implemented'.format(type))
        interpolatesv.requires_grad_(True)
        disc_interpolates = netD(interpolatesv)
        gradients = torch.autograd.grad(outputs=disc_interpolates, inputs=interpolatesv,
                                        grad_outputs=torch.ones(disc_interpolates.size()).to(device),
                                        create_graph=True, retain_graph=True, only_inputs=True)
        gradients = gradients[0].view(real_data.size(0), -1)  # flat the data
        gradient_penalty = (((gradients + 1e-16).norm(2, dim=1) - constant) ** 2).mean() * lambda_gp  # added eps
        return gradient_penalty, gradients
    else:
        return 0.0, None


class ResnetGenerator(nn.Module):
    """Resnet-based generator that consists of Resnet blocks between a few downsampling/upsampling operations.
    We adapt Torch code and idea from Justin Johnson's neural style transfer project(https://github.com/jcjohnson/fast-neural-style)
    """

    def __init__(self, input_nc, output_nc, ngf=64, n_downsampling=2, norm_layer=nn.BatchNorm2d, use_dropout=True,
                 n_blocks=6,
                 padding_type='zero'):
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
        assert (n_blocks >= 0)
        super(ResnetGenerator, self).__init__()
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


class ResnetBlock(nn.Module):
    """Define a Resnet block"""

    def __init__(self, dim, padding_type, norm_layer, use_dropout, use_bias):
        """Initialize the Resnet block

        A resnet block is a conv block with skip connections
        We construct a conv block with build_conv_block function,
        and implement skip connections in <forward> function.
        Original Resnet paper: https://arxiv.org/pdf/1512.03385.pdf
        """
        super(ResnetBlock, self).__init__()
        self.conv_block = self.build_conv_block(dim, padding_type, norm_layer, use_dropout, use_bias)

    def build_conv_block(self, dim, padding_type, norm_layer, use_dropout, use_bias):
        """Construct a convolutional block.

        Parameters:
            dim (int)           -- the number of channels in the conv layer.
            padding_type (str)  -- the name of padding layer: reflect | replicate | zero
            norm_layer          -- normalization layer
            use_dropout (bool)  -- if use dropout layers.
            use_bias (bool)     -- if the conv layer uses bias or not

        Returns a conv block (with a conv layer, a normalization layer, and a non-linearity layer (ReLU))
        """
        conv_block = []
        p = 0
        if padding_type == 'reflect':
            conv_block += [nn.ReflectionPad2d(1)]
        elif padding_type == 'replicate':
            conv_block += [nn.ReplicationPad2d(1)]
        elif padding_type == 'zero':
            p = 1
        else:
            raise NotImplementedError('padding [%s] is not implemented' % padding_type)

        conv_block += [nn.Conv2d(dim, dim, kernel_size=3, padding=p, bias=use_bias), norm_layer(dim), nn.ReLU(True)]
        if use_dropout:
            conv_block += [nn.Dropout(0.25)]  # original dropout = 0.5

        p = 0
        if padding_type == 'reflect':
            conv_block += [nn.ReflectionPad2d(1)]
        elif padding_type == 'replicate':
            conv_block += [nn.ReplicationPad2d(1)]
        elif padding_type == 'zero':
            p = 1
        else:
            raise NotImplementedError('padding [%s] is not implemented' % padding_type)
        conv_block += [nn.Conv2d(dim, dim, kernel_size=3, padding=p, bias=use_bias), norm_layer(dim)]

        return nn.Sequential(*conv_block)

    def forward(self, x):
        """Forward function (with skip connections)"""
        out = x + self.conv_block(x)  # add skip connections
        return out


class UnetGenerator(nn.Module):
    """Create a Unet-based generator"""

    def __init__(self, input_nc, output_nc, num_downs, ngf=64, norm_layer=nn.BatchNorm2d, use_dropout=False,
                 attention=False):
        """Construct a Unet generator
        Parameters:
            input_nc (int)  -- the number of channels in input images
            output_nc (int) -- the number of channels in output images
            num_downs (int) -- the number of downsamplings in UNet. For example, # if |num_downs| == 7,
                                image of size 128x128 will become of size 1x1 # at the bottleneck
            ngf (int)       -- the number of filters in the last conv layer
            norm_layer      -- normalization layer
        """
        super(UnetGenerator, self).__init__()
        # construct unet structure
        unet_block = UnetSkipConnectionBlock(ngf * 8, ngf * 8, input_nc=None, submodule=None, norm_layer=norm_layer,
                                             innermost=True, use_dropout=True)  # add the innermost layer
        for i in range(num_downs - 5):  # add intermediate layers with ngf * 8 filters
            unet_block = UnetSkipConnectionBlock(ngf * 8, ngf * 8, input_nc=None, submodule=unet_block,
                                                 norm_layer=norm_layer, use_dropout=use_dropout)
        # gradually reduce the number of filters from ngf * 8 to ngf
        unet_block = UnetSkipConnectionBlock(ngf * 4, ngf * 8, input_nc=None, submodule=unet_block,
                                             norm_layer=norm_layer)
        unet_block = UnetSkipConnectionBlock(ngf * 2, ngf * 4, input_nc=None, submodule=unet_block,
                                             norm_layer=norm_layer)
        unet_block = UnetSkipConnectionBlock(ngf, ngf * 2, input_nc=None, submodule=unet_block, norm_layer=norm_layer)
        self.model = UnetSkipConnectionBlock(output_nc, ngf, input_nc=input_nc, submodule=unet_block, outermost=True,
                                             norm_layer=norm_layer, attention=attention)  # add the outermost layer

    def forward(self, input):
        """Standard forward"""
        return self.model(input)


class UnetSkipConnectionBlock(nn.Module):
    """Defines the Unet submodule with skip connection.
         -------------------identity----------------------
        |-- downsampling -- |submodule| -- upsampling --|
    """

    def __init__(self, outer_nc, inner_nc, input_nc=None,
                 submodule=None, outermost=False, innermost=False, norm_layer=nn.BatchNorm2d, use_dropout=False,
                 attention=False):
        """Construct a Unet submodule with skip connections.

        Parameters:
            outer_nc (int) -- the number of filters in the outer conv layer
            inner_nc (int) -- the number of filters in the inner conv layer
            input_nc (int) -- the number of channels in input images/features
            submodule (UnetSkipConnectionBlock) -- previously defined submodules
            outermost (bool)    -- if this module is the outermost module
            innermost (bool)    -- if this module is the innermost module
            norm_layer          -- normalization layer
            use_dropout (bool)  -- if use dropout layers.
        """
        super(UnetSkipConnectionBlock, self).__init__()
        self.outermost = outermost
        if type(norm_layer) == functools.partial:
            use_bias = norm_layer.func == nn.InstanceNorm2d
        else:
            use_bias = norm_layer == nn.InstanceNorm2d
        if input_nc is None:
            input_nc = outer_nc
        downconv = nn.Conv2d(input_nc, inner_nc, kernel_size=4,
                             stride=2, padding=1, bias=use_bias)
        downrelu = nn.LeakyReLU(0.2, True)
        downnorm = norm_layer(inner_nc)
        uprelu = nn.ReLU(True)
        upnorm = norm_layer(outer_nc)

        if outermost:
            upconv = nn.ConvTranspose2d(inner_nc * 2, outer_nc,
                                        kernel_size=4, stride=2,
                                        padding=1)
            if attention:
                down_sa = Self_Attn(input_nc)
                up_sa = Self_Attn(inner_nc * 2)
                down = [down_sa, downconv]
                up = [uprelu, up_sa, upconv, nn.Tanh()]
            else:
                down = [downconv]
                up = [uprelu, upconv, nn.Tanh()]
            model = down + [submodule] + up
        elif innermost:
            upconv = nn.ConvTranspose2d(inner_nc, outer_nc,
                                        kernel_size=4, stride=2,
                                        padding=1, bias=use_bias)
            down = [downrelu, downconv]
            up = [uprelu, upconv, upnorm]
            model = down + up
        else:
            upconv = nn.ConvTranspose2d(inner_nc * 2, outer_nc,
                                        kernel_size=4, stride=2,
                                        padding=1, bias=use_bias)

            # Self-Attn
            if attention:
                down_sa = Self_Attn(input_nc)
                up_sa = Self_Attn(inner_nc * 2)
                down = [downrelu, down_sa, downconv, downnorm]
                up = [uprelu, up_sa, upconv, upnorm]
            else:
                down = [downrelu, downconv, downnorm]
                up = [uprelu, upconv, upnorm]

            # without SA
            # down = [downrelu, downconv, downnorm]
            # up = [uprelu, upconv, upnorm]

            if use_dropout:
                model = down + [submodule] + up + [nn.Dropout(0.25)]
            else:
                model = down + [submodule] + up

        self.model = nn.Sequential(*model)

    def forward(self, x):
        if self.outermost:
            return self.model(x)
        else:  # add skip connections
            return torch.cat([x, self.model(x)], 1)


class SimpleEnhancer(nn.Module):
    def __init__(self, input_nc, output_nc, ngf=32, n_downsample_global=2, n_blocks_global=9,
                 n_local_enhancers=1, n_blocks_local=9, norm_layer=nn.BatchNorm2d, padding_type='zero'):
        super(SimpleEnhancer, self).__init__()
        self.n_local_enhancers = n_local_enhancers

        ###### global generator model #####
        ngf_global = ngf * (2 ** n_local_enhancers)
        model_global = ResnetGenerator(input_nc, output_nc, ngf_global, n_downsample_global, n_blocks=n_blocks_global,
                                       norm_layer=norm_layer).model
        model_global = [model_global[i] for i in
                        range(len(model_global) - 3)]  # get rid of final convolution layers
        self.model = nn.Sequential(*model_global)

        ###### local enhancer layers #####
        for n in range(1, n_local_enhancers + 1):
            ### downsample
            ngf_global = ngf * (2 ** (n_local_enhancers - n))
            model_downsample = [nn.ZeroPad2d(1), nn.Conv2d(input_nc, ngf_global, kernel_size=3, padding=0),
                                norm_layer(ngf_global), nn.ReLU(True),
                                nn.Conv2d(ngf_global, ngf_global * 2, kernel_size=3, stride=2, padding=1),
                                norm_layer(ngf_global * 2), nn.ReLU(True)
                                ]
            ### residual blocks
            model_upsample = []
            for i in range(n_blocks_local):
                model_upsample += [ResnetBlock(ngf_global * 2, padding_type=padding_type, norm_layer=norm_layer,
                                               use_bias=True, use_dropout=True)]

            ### upsample
            model_upsample += [
                nn.ConvTranspose2d(ngf_global * 2, ngf_global, kernel_size=3, stride=2, padding=1, output_padding=1),
                norm_layer(ngf_global), nn.ReLU(True)]

            ### final convolution
            if n == n_local_enhancers:
                model_upsample += [nn.ZeroPad2d(1), nn.Conv2d(ngf, output_nc, kernel_size=3, padding=0),
                                   nn.Tanh()]

            setattr(self, 'model' + str(n) + '_1', nn.Sequential(*model_downsample))
            setattr(self, 'model' + str(n) + '_2', nn.Sequential(*model_upsample))

        self.downsample = nn.Sequential(
            # nn.AvgPool2d(3, stride=2, padding=[1, 1], count_include_pad=False),
            nn.MaxPool2d(2)
        )

    def forward(self, input):
        ### create input pyramid
        input_downsampled = [input]
        for i in range(self.n_local_enhancers):
            input_downsampled.append(self.downsample(input_downsampled[-1]))

        ### output at coarest level
        output_prev = self.model(input_downsampled[-1])
        ### build up one layer at a time
        for n_local_enhancers in range(1, self.n_local_enhancers + 1):
            model_downsample = getattr(self, 'model' + str(n_local_enhancers) + '_1')
            model_upsample = getattr(self, 'model' + str(n_local_enhancers) + '_2')
            input_i = input_downsampled[self.n_local_enhancers - n_local_enhancers]
            output_prev = model_upsample(model_downsample(input_i) + output_prev)
        return output_prev


class Filling(nn.Module):
    def __init__(self, input_nc, output_nc, ngf=32, n_blocks=9, norm_layer=nn.BatchNorm2d, padding_type='zero'):
        super(Filling, self).__init__()

        ##### filling layers #####
        ### downsample
        model_downsample = [nn.ZeroPad2d(1), nn.Conv2d(input_nc, ngf, kernel_size=3, padding=0),
                            norm_layer(ngf), nn.ReLU(True),
                            nn.Conv2d(ngf, ngf * 2, kernel_size=3, stride=2, padding=1),
                            norm_layer(ngf * 2), nn.ReLU(True)
                            ]
        ### residual blocks
        model_upsample = []
        for i in range(n_blocks):
            model_upsample += [ResnetBlock(ngf * 2, padding_type=padding_type, norm_layer=norm_layer,
                                           use_bias=True, use_dropout=True)]

        ### upsample
        model_upsample += [
            nn.ConvTranspose2d(ngf * 2, ngf, kernel_size=3, stride=2, padding=1, output_padding=1),
            norm_layer(ngf), nn.ReLU(True)]

        ### final convolution
        model_upsample += [nn.ZeroPad2d(1), nn.Conv2d(ngf, output_nc, kernel_size=3, padding=0),
                           nn.Tanh()]

        setattr(self, 'model' + str(1) + '_1', nn.Sequential(*model_downsample))
        setattr(self, 'model' + str(1) + '_2', nn.Sequential(*model_upsample))

    def forward(self, input):
        model_downsample = getattr(self, 'model' + str(1) + '_1')
        model_upsample = getattr(self, 'model' + str(1) + '_2')
        output = model_upsample(model_downsample(input))
        return output


class Discriminator(nn.Module):
    def __init__(self, input_nc, ndf=64, n_layers=3, norm_layer=nn.BatchNorm2d, use_dropout=False):
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


class Critic(nn.Module):
    def __init__(self, input_nc, ndf=64, n_layers=3, norm_layer=nn.BatchNorm2d, use_dropout=False):
        """Construct a PatchGAN discriminator

        Parameters:
            input_nc (int)  -- the number of channels in input images
            ndf (int)       -- the number of filters in the last conv layer
            n_layers (int)  -- the number of conv layers in the discriminator
            norm_layer      -- normalization layer
        """
        super(Critic, self).__init__()
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

        self.constraint_net = nn.Sequential(
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
        correctness = self.constraint_net(preceded_fm)
        return correctness


class Multi_Discriminator(nn.Module):
    def __init__(self, input_nc, ndf=64, n_layers=3, norm_layer=nn.BatchNorm2d, use_dropout=False):
        """Construct a Bi-discriminator.

        Parameters:
            input_nc (int)  -- the number of channels in input images
            ndf (int)       -- the number of filters in the last conv layer
            n_layers (int)  -- the number of conv layers in the discriminator
            norm_layer      -- normalization layer
        """
        super(Multi_Discriminator, self).__init__()
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

        self.constraint_net = nn.Sequential(
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
        correctness = self.constraint_net(preceded_fm)
        return out, correctness


class Self_Attn(nn.Module):
    """ Self attention Layer"""

    def __init__(self, in_dim):
        super(Self_Attn, self).__init__()
        self.channel_in = in_dim

        self.query_conv = nn.Conv2d(in_channels=in_dim, out_channels=8, kernel_size=1)
        self.key_conv = nn.Conv2d(in_channels=in_dim, out_channels=8, kernel_size=1)
        self.value_conv = nn.Conv2d(in_channels=in_dim, out_channels=in_dim, kernel_size=1)
        self.gamma = nn.Parameter(torch.zeros(1))

        self.softmax = nn.Softmax(dim=-1)  #

    def forward(self, x):
        """
            inputs :
                x : input feature maps( B X C X W X H)
            returns :
                out : self attention value + input feature
                attention: B X N X N (N is Width*Height)
        """
        m_batchsize, C, width, height = x.size()
        proj_query = self.query_conv(x).view(m_batchsize, -1, width * height).permute(0, 2, 1)  # B X CX(N)
        proj_key = self.key_conv(x).view(m_batchsize, -1, width * height)  # B X C x (*W*H)
        energy = torch.bmm(proj_query, proj_key)  # transpose check
        attention = self.softmax(energy)  # BX (N) X (N)
        proj_value = self.value_conv(x).view(m_batchsize, -1, width * height)  # B X C X N

        out = torch.bmm(proj_value, attention.permute(0, 2, 1))
        out = out.view(m_batchsize, C, width, height)

        out = self.gamma * out + x
        return out


class CVAE(nn.Module):
    """ Construct a generator which aims at reproducing the model proposed by the 2-page paper in DAC'20. """
    def __init__(self, input_nc, output_nc, latent_dim=128, ngf=32, norm_layer=nn.BatchNorm2d,
                 use_dropout=True):
        super(CVAE, self).__init__()

        self.encoder = nn.Sequential(
            nn.Conv2d(input_nc, ngf, kernel_size=3, stride=1, padding=1, bias=True),
            norm_layer(ngf), nn.ReLU(True),
            nn.AvgPool2d(2),
            nn.Conv2d(ngf, ngf * 2, kernel_size=3, stride=1, padding=1, bias=True),
            norm_layer(ngf * 2), nn.ReLU(True),
            nn.AvgPool2d(2),
            nn.Conv2d(ngf * 2, ngf * 4, kernel_size=3, stride=1, padding=1, bias=True),
            norm_layer(ngf * 4), nn.ReLU(True),
            nn.AvgPool2d(2),
            nn.Conv2d(ngf * 4, ngf * 8, kernel_size=3, stride=1, padding=1, bias=True),
            norm_layer(ngf * 8), nn.ReLU(True),
            nn.AvgPool2d(2),
            nn.Conv2d(ngf * 8, ngf * 8, kernel_size=4, stride=1, padding=0, bias=True),
            norm_layer(ngf * 8), nn.ReLU(True)
        )

        self.dense = nn.Sequential(
            nn.Linear(ngf * 8, ngf * 16),
            nn.Linear(ngf * 16, ngf * 8)
        )

        self.decoder = nn.Sequential(
            nn.ConvTranspose2d(ngf * 8, ngf * 8, kernel_size=4, stride=1, padding=0, bias=True),
            norm_layer(ngf * 8), nn.ReLU(True),
            nn.ConvTranspose2d(ngf * 8, ngf * 4, kernel_size=4, stride=2, padding=1, bias=True),
            norm_layer(ngf * 4), nn.ReLU(True),
            nn.ConvTranspose2d(ngf * 4, ngf * 2, kernel_size=4, stride=2, padding=1, bias=True),
            norm_layer(ngf * 2), nn.ReLU(True),
            nn.ConvTranspose2d(ngf * 2, ngf, kernel_size=4, stride=2, padding=1, bias=True),
            norm_layer(ngf), nn.ReLU(True),
            nn.ConvTranspose2d(ngf, output_nc, kernel_size=4, stride=2, padding=1, bias=True),
            nn.Tanh()
        )

    def forward(self, x):
        """Standard forward"""
        x = self.encoder(x).view(x.size(0), -1)
        x = self.dense(x).unsqueeze(-1).unsqueeze(-1)
        out = self.decoder(x)
        return out
