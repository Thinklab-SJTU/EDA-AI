import torch
import numpy as np
from . import networks
from .base_model import BaseModel
from . import L1NormLoss
from utils import utils
from .loss import *


class CVAEModel(BaseModel):
    """ This class implements the CVAE model, for learning a mapping from input images to output images given paired
    data.

    The model training requires '--dataset_mode aligned' dataset.
    By default, it uses a '--netG cvae' generator.
    """

    @staticmethod
    def modify_commandline_options(parser, is_train=True):
        """Add new dataset-specific options, and rewrite default values for existing options.

        Parameters:
            parser          -- original option parser
            is_train (bool) -- whether training phase or test phase. You can use this flag to add training-specific or test-specific options.

        Returns:
            the modified parser.

        For pix2pix, we do not use image buffer
        The training objective is: GAN Loss + lambda_L1 * ||G(A)-B||_1
        By default, we use vanilla GAN loss, UNet with batchnorm, and aligned datasets.
        """
        if is_train:
            parser.add_argument('--lambda_focal', type=float, default=100, help='weight for focal loss')
            parser.add_argument('--lambda_norm', type=float, default=100, help='weight for norm loss')
            parser.add_argument('--mix', type=float, default=0.1,
                                help='weight for nir select optimal layer L1 loss')

        return parser

    def __init__(self, opt):
        """Initialize the class.

        Parameters:
            opt (Option class)-- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        BaseModel.__init__(self, opt)
        # specify the training losses you want to print out. The training/test scripts will call
        # <BaseModel.get_current_losses>
        self.loss_names = ['G_norm']
        # specify the images you want to save/display. The training/test scripts will call
        # <BaseModel.get_current_visuals>
        self.visual_names = ['real_A', 'fake_B', 'real_B']
        # specify the models you want to save to the disk. The training/test scripts will call
        # <BaseModel.save_networks> and <BaseModel.load_networks>

        if self.isTrain:
            self.model_names = ['G']
        else:  # during test time, only load G
            self.model_names = ['G']
        # define networks (both generator and discriminator)
        self.netG = networks.define_G(opt.input_nc, opt.output_nc, opt.ngf, opt.netG, opt.norm,
                                      not opt.no_dropout, opt.init_type, opt.init_gain, self.gpu_ids)

        if self.isTrain:
            # define loss functions
            self.criterionGAN = networks.GANLoss(opt.gan_mode).to(self.device)
            self.criterionMSE = torch.nn.MSELoss()
            # initialize optimizers; schedulers will be automatically created by function <BaseModel.setup>.
            self.optimizer_G = torch.optim.Adam(self.netG.parameters(), lr=opt.lr, betas=(opt.beta1, 0.999),
                                                weight_decay=0.01)
            self.optimizers.append(self.optimizer_G)

    def set_input(self, input):
        """Unpack input data from the dataloader and perform necessary pre-processing steps.

        Parameters:
            input (dict): include the data itself and its metadata information.
        """
        self.real_A = input['A'].to(self.device)
        self.real_B = input['B'].to(self.device)
        self.real_label = (self.real_B + 1.) / 2.
        # self.image_paths = input['A_paths']

    def forward(self):
        """Run forward pass; called by both functions <optimize_parameters> and <test>."""
        self.fake_B = self.netG(self.real_A)  # G(A)

        # record the legally generated images in the current epoch
        self.current_correctness = utils.constraint_label(self.real_A.cpu(),
                                                          self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        # compute the total effective length of labeled and generated routes
        self.realB_length = (self.real_B / 2 + 0.5).view(self.real_B.size(0), -1).sum(1).unsqueeze(-1)
        self.fake_length = torch.round(self.fake_B / 2 + 0.5).view(self.fake_B.size(0), -1).sum(1).unsqueeze(-1)

        # aggregate the generated images which satisfy the rules
        self.correctness += self.current_correctness.sum().item()
        self.aligned_real_twl += (self.current_correctness * self.realB_length).sum().item()
        self.aligned_fake_twl += (self.current_correctness * self.fake_length).sum().item()
        self.dist = (self.fake_length - self.realB_length).sum().item()
        self.current_correctness = self.current_correctness.unsqueeze(-1).unsqueeze(-1)

    def backward_G(self):
        """Calculate the weighted MSE loss for the generator"""
        # Second, G(A) = B
        self.loss_G_norm = self.criterionMSE(self.fake_B, self.real_B)

        # combine loss and calculate gradients
        self.loss_G = self.loss_G_norm * (1 + 1e-3 * (100 * np.sign(self.dist - 1) + 1) *
                                          self.dist)
        self.loss_G.backward()

    def optimize_parameters(self):
        self.forward()  # compute fake images: G(A)
        # update G
        self.optimizer_G.zero_grad()  # set G's gradients to zero
        self.backward_G()  # calculate graidents for G
        self.optimizer_G.step()  # udpate G's weights
