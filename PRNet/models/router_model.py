import numpy as np
import torch
from . import networks
from .base_model import BaseModel
from . import L1NormLoss
from utils import utils
from .loss import *


class RouterModel(BaseModel):
    """ This class implements the pix2pix model, for learning a mapping from input images to output images given paired
    data.

    The model training requires '--dataset_mode aligned' dataset.
    By default, it uses a '--netG unet_512' U-Net generator,
    a '--netD basic' discriminator (PatchGAN),
    and a '--gan_mode' vanilla GAN loss (the cross-entropy objective used in the orignal GAN paper).

    pix2pix paper: https://arxiv.org/pdf/1611.07004.pdf
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
        # changing the default values to match the pix2pix paper (https://phillipi.github.io/pix2pix/)
        # parser.set_defaults(norm='batch', netG='unet_512', dataset_mode='aligned')
        if is_train:
            # parser.set_defaults(pool_size=0, gan_mode='vanilla')
            parser.add_argument('--lambda_focal', type=float, default=100, help='weight for focal loss')
            parser.add_argument('--lambda_norm', type=float, default=100, help='weight for norm loss')
            parser.add_argument('--lambda_L1', type=float, default=100.0, help='weight for L1 loss')
            parser.add_argument('--mix', type=float, default=0.1,
                                help='weight for nir select optimal layer L1 loss')

        return parser

    def __init__(self, opt):
        """Initialize the pix2pix class.

        Parameters:
            opt (Option class)-- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        BaseModel.__init__(self, opt)
        # specify the training losses you want to print out. The training/test scripts will call
        # <BaseModel.get_current_losses>
        # self.loss_names = ['G_GAN', 'G_norm', 'G_focal', 'D_real', 'D_fake', 'D_constraint', 'D_jam', 'D_length']
        self.loss_names = ['G_GAN', 'G_L1', 'D_real', 'D_fake']
        # specify the images you want to save/display. The training/test scripts will call
        # <BaseModel.get_current_visuals>
        self.visual_names = ['real_A', 'fake_B', 'real_B']
        # specify the models you want to save to the disk. The training/test scripts will call
        # <BaseModel.save_networks> and <BaseModel.load_networks>

        if self.isTrain:
            self.model_names = ['G', 'D']
        else:  # during test time, only load G
            self.model_names = ['G']
        # define networks (both generator and discriminator)
        self.netG = networks.define_G(opt.input_nc, opt.output_nc, opt.ngf, opt.netG, opt.norm,
                                      not opt.no_dropout, opt.init_type, opt.init_gain, self.gpu_ids)

        if self.isTrain:  # define a discriminator; conditional GANs need to take both input and output images;
            # Therefore, #channels for D is input_nc + output_nc
            self.netD = networks.define_D(opt.input_nc + opt.output_nc, opt.ndf, opt.netD,
                                          opt.n_layers_D, opt.norm, opt.init_type, opt.init_gain, self.gpu_ids)

        if self.isTrain:
            # define loss functions
            self.criterionGAN = networks.GANLoss(opt.gan_mode).to(self.device)
            self.criterionMSE = torch.nn.MSELoss()
            self.criterionBCE = torch.nn.BCELoss()
            self.focalloss = FocalLoss()
            self.normloss = NormLoss(2)
            self.criterionL1 = torch.nn.L1Loss()
            self.criterionL1_M = L1NormLoss.L1_norm()
            # initialize optimizers; schedulers will be automatically created by function <BaseModel.setup>.
            self.optimizer_G = torch.optim.Adam(self.netG.parameters(), lr=opt.lr, betas=(opt.beta1, 0.999),
                                                weight_decay=0.01)
            self.optimizer_D = torch.optim.Adam(self.netD.parameters(), lr=opt.lr, betas=(opt.beta1, 0.999),
                                                weight_decay=0.01)
            self.optimizers.append(self.optimizer_G)
            self.optimizers.append(self.optimizer_D)

    def set_input(self, input):
        """Unpack input data from the dataloader and perform necessary pre-processing steps.

        Parameters:
            input (dict): include the data itself and its metadata information.

        The option 'direction' can be used to swap images in domain A and domain B.
        """
        self.real_A = input['A'].to(self.device)
        self.real_B = input['B'].to(self.device)
        self.image_paths = input['A_paths']

    def forward(self):
        """Run forward pass; called by both functions <optimize_parameters> and <test>."""
        self.fake_B = self.netG(self.real_A)  # G(A)
        self.current_correctness = utils.constraint_label(self.real_A.cpu(),
                                                          self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        realB_length = utils.length_label(self.real_B.cpu()).unsqueeze(-1).to(self.device)
        fake_length = utils.length_label(self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        self.correctness += self.current_correctness.sum().item()
        self.aligned_real_twl += (self.current_correctness * realB_length).sum().item()
        self.aligned_fake_twl += (self.current_correctness * fake_length).sum().item()
        self.dist = (fake_length - realB_length).sum().item()

        # self.fake_length = utils.length_label(self.fake_B.cpu()).unsqueeze(-1).to(self.device)

    def backward_D(self):
        """Calculate GAN loss for the discriminator"""
        # Fake; stop backprop to the generator by detaching fake_B
        # fake_AB = torch.cat((self.real_A, self.fake_B), 1)  # we use conditional GANs; we need to feed both input and
        # # output to the discriminator
        # pred_fake, pred_constraint, pred_jam, pred_length = self.netD(fake_AB.detach())
        # self.loss_D_fake = self.criterionGAN(pred_fake, False)
        #
        # constraint_real = utils.constraint_label(self.real_A.cpu(), self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        # self.connectivity += constraint_real.sum().item()
        # self.loss_D_constraint = self.criterionBCE(pred_constraint, constraint_real)
        #
        # jam_real = utils.jam_label(self.real_A.cpu(), self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        # self.congestion += jam_real.sum().item()
        # self.loss_D_jam = self.criterionMSE(pred_jam, jam_real)
        #
        # self.length_real = utils.length_label(self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        # self.wirelength += self.length_real.sum().item()
        # self.loss_D_length = self.criterionMSE(pred_length, self.length_real)
        #
        # # Real
        # real_AB = torch.cat((self.real_A, self.real_B), 1)
        # pred_real, _, _, _ = self.netD(real_AB)
        # self.loss_D_real = self.criterionGAN(pred_real, True)
        # # self.loss_D_constraint += self.criterionGAN(pred_constraint, True)
        # # combine loss and calculate gradients
        # self.loss_D = self.loss_D_fake + self.loss_D_real + self.loss_D_constraint + self.loss_D_jam + self.loss_D_length
        # self.loss_D.backward(retain_graph=True)

        fake_AB = torch.cat((self.real_A, self.fake_B), 1)  # we use conditional GANs; we need to feed both input and
        # output to the discriminator
        pred_fake = self.netD(fake_AB.detach())
        self.loss_D_fake = self.criterionGAN(pred_fake, False)
        # Real
        real_AB = torch.cat((self.real_A, self.real_B), 1)
        pred_real = self.netD(real_AB)
        self.loss_D_real = self.criterionGAN(pred_real, True)
        # combine loss and calculate gradients
        self.loss_D = (self.loss_D_fake + self.loss_D_real) * 0.5
        self.loss_D.backward()

        # print(pred_fake.mean(), pred_real.mean())

    def backward_G(self):
        """Calculate GAN and L1 loss for the generator"""
        # First, G(A) should fake the discriminator
        # fake_AB = torch.cat((self.real_A, self.fake_B), 1)
        # pred_fake, _, _, _ = self.netD(fake_AB)
        # realB_length = utils.length_label(self.real_B.cpu()).unsqueeze(-1).to(self.device)
        # self.loss_G_GAN = self.criterionGAN(pred_fake, True)
        # # Second, G(A) = B
        # self.loss_G_norm = self.normloss(self.fake_B, self.real_B) * self.opt.lambda_norm
        # self.loss_G_focal = self.focalloss((self.fake_B + 1) / 2, self.real_label) * self.opt.lambda_focal
        # # self.loss_M_L1 = self.criterionL1_M(self.M_L1)
        # # print(self.loss_M_L1)
        # # combine loss and calculate gradients
        # coef = torch.abs(self.length_real - realB_length).sum().item()
        # self.loss_G = (self.loss_G_norm + self.loss_G_focal + self.loss_G_GAN) * (1 + coef)  # + self.loss_M_L1*0.1
        # # self.loss_G = (self.loss_G_norm + self.loss_G_focal + self.loss_G_GAN)
        # # print(self.loss_G)
        # self.loss_G.backward()

        """Calculate GAN and L1 loss for the generator"""
        # First, G(A) should fake the discriminator
        fake_AB = torch.cat((self.real_A, self.fake_B), 1)
        pred_fake = self.netD(fake_AB)
        self.loss_G_GAN = self.criterionGAN(pred_fake, True)
        # Second, G(A) = B

        # realB_length = utils.length_label(self.real_B.cpu()).unsqueeze(-1).to(self.device)
        # fake_length = utils.length_label(self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        self.loss_G_L1 = self.criterionMSE(self.fake_B, self.real_B) * (1 + 1e-3 * (100 * np.sign(self.dist - 1) + 1) *
                                                                        self.dist)
        # self.loss_G_L1 = self.criterionMSE(self.fake_B, self.real_B) * self.opt.lambda_L1

        # correctness = utils.constraint_label(self.real_A.cpu(), self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        # self.correctness += correctness.sum().item()
        # self.aligned_real_twl += (correctness * realB_length).sum().item()
        # self.aligned_fake_twl += (correctness * fake_length).sum().item()


        # combine loss and calculate gradients
        self.loss_G = self.loss_G_GAN + self.loss_G_L1
        self.loss_G.backward()

    def optimize_parameters(self):
        self.forward()  # compute fake images: G(A)
        # update D
        self.set_requires_grad(self.netD, True)  # enable backprop for D
        self.optimizer_D.zero_grad()  # set D's gradients to zero
        self.backward_D()  # calculate gradients for D
        self.optimizer_D.step()  # update D's weights
        # update G
        self.set_requires_grad(self.netD, False)  # D requires no gradients when optimizing G
        self.optimizer_G.zero_grad()  # set G's gradients to zero
        self.backward_G()  # calculate graidents for G
        self.optimizer_G.step()  # udpate G's weights'




