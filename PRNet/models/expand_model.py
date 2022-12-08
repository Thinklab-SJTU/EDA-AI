import torch
import numpy as np
from . import networks
from .base_model import BaseModel
from . import L1NormLoss
from utils import utils
from .loss import *


class ExpandModel(BaseModel):
    """ This class implements the Expanded model, for learning a mapping from input images to output images given paired
    data.

    The model training requires '--dataset_mode aligned' dataset.
    By default, it uses a '--netG g_128_f' generator,
    a '--netD basic' discriminator,
    and a '--gan_mode' vanilla GAN loss (the cross-entropy objective used in the orignal GAN paper).
    """

    @staticmethod
    def modify_commandline_options(parser, is_train=True):
        """Add new dataset-specific options, and rewrite default values for existing options.

        Parameters:
            parser          -- original option parser
            is_train (bool) -- whether training phase or test phase. You can use this flag to add training-specific or test-specific options.

        Returns:
            the modified parser.
        """
        if is_train:
            parser.add_argument('--lambda_focal', type=float, default=10, help='weight for focal loss')
            parser.add_argument('--lambda_norm', type=float, default=10, help='weight for norm loss')
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
        self.loss_names = ['G_GAN', 'G_norm', 'G_focal', 'D_real', 'D_fake', 'C_fake', 'C_real']
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
        self.netG = networks.define_G(opt.input_nc, opt.output_nc, opt.ngf, 'g_128_f', opt.norm,
                                      not opt.no_dropout, opt.init_type, opt.init_gain, self.gpu_ids)

        if self.isTrain:  # define a discriminator; conditional GANs need to take both input and output images;
            # Therefore, #channels for D is input_nc + output_nc
            self.netD = networks.define_D(opt.input_nc + opt.output_nc, opt.ndf, opt.netD, opt.n_layers_D, opt.norm,
                                          opt.init_type, opt.init_gain, self.gpu_ids)

        if self.isTrain:
            # define loss functions
            self.criterionGAN = networks.GANLoss(opt.gan_mode).to(self.device)
            self.criterionMSE = torch.nn.MSELoss()
            self.criterionBCE = torch.nn.BCEWithLogitsLoss()
            self.focalloss = FocalLoss()
            self.normloss = NormLoss(2)
            self.criterionL1 = torch.nn.L1Loss()
            # initialize optimizers; schedulers will be automatically created by function <BaseModel.setup>.
            self.optimizer_G = torch.optim.Adam(self.netG.parameters(), lr=opt.lr, betas=(opt.beta1, 0.999),
                                                weight_decay=0.01)
            self.optimizer_D = torch.optim.Adam(self.netD.parameters(), lr=opt.lr, betas=(opt.beta1, 0.999),
                                                weight_decay=0.01)
            self.optimizers.append(self.optimizer_G)
            self.optimizers.append(self.optimizer_D)

    def initialize(self, opt):
        if self.isTrain:
            # load partly the trained basic generator network as the guiding network, and the pre-trained filling network
            self.load_subNetworks(self.netG, 'G', epoch, pretrained_path)
            self.load_subNetworks(self.netG, 'G', epoch, pretrained_path)
            self.load_subNetworks(self.netD, 'D', epoch, pretrained_path)
            self.schedulers = [networks.get_scheduler(optimizer, opt) for optimizer in self.optimizers]
        if not self.isTrain or opt.continue_train:
            load_suffix = 'iter_%d' % opt.load_iter if opt.load_iter > 0 else opt.epoch
            self.load_networks(load_suffix)
        self.print_networks(opt.verbose)

    def set_input(self, input):
        """Unpack input data from the dataloader and perform necessary pre-processing steps.

        Parameters:
            input (dict): include the data itself and its metadata information.
        """
        self.real_A = input['A'].to(self.device)
        self.real_B = input['B'].to(self.device)
        self.real_label = (self.real_B + 1.) / 2.

    def forward(self):
        """Run forward pass; called by both functions <optimize_parameters> and <test>."""
        self.fake_B = self.netG(self.real_A)  # G(A)

        self.current_correctness = utils.constraint_label(self.real_A.cpu(),
                                                          self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        self.realB_length = (self.real_B / 2 + 0.5).view(self.real_B.size(0), -1).sum(1).unsqueeze(-1)
        self.fake_length = torch.round(self.fake_B / 2 + 0.5).view(self.fake_B.size(0), -1).sum(1).unsqueeze(-1)

        self.correctness += self.current_correctness.sum().item()
        self.aligned_real_twl += (self.current_correctness * self.realB_length).sum().item()
        self.aligned_fake_twl += (self.current_correctness * self.fake_length).sum().item()
        self.dist = (self.fake_length - self.realB_length).sum().item()
        self.current_correctness = self.current_correctness.unsqueeze(-1).unsqueeze(-1)

        self.hpwl = utils.hpwl(self.real_A.cpu()).unsqueeze(-1).to(self.device)

    def backward_D(self):
        """Calculate GAN loss for the discriminator"""
        # Fake: stop backprop to the generator by detaching fake_B
        fake_AB = torch.cat((self.real_A, self.fake_B), 1)  # we use conditional GANs; we need to feed both input and
        # output to the discriminator
        pred_fake, pred_fc = self.netD(fake_AB.detach())
        self.loss_D_fake = self.criterionGAN(pred_fake, False)
        self.loss_C_fake = self.criterionBCE(pred_fc, self.current_correctness.expand_as(pred_fc))

        # Real
        real_AB = torch.cat((self.real_A, self.real_B), 1)
        pred_real, pred_rc = self.netD(real_AB)
        self.loss_D_real = self.criterionGAN(pred_real, True)
        self.loss_C_real = self.criterionGAN(pred_rc, True)

        # combine loss and calculate gradients
        self.loss_D = (self.loss_D_fake + self.loss_D_real) * 0.5
        self.loss_C = (self.loss_C_fake + self.loss_C_real) * 0.5
        self.loss_D.backward(retain_graph=True)
        self.loss_C.backward()

        self.D_real_loss += self.loss_D_real.sum().item()
        self.D_fake_loss += self.loss_D_fake.sum().item()
        self.D_C += self.loss_C.sum().item()

    def backward_G(self):
        """Calculate GAN and Norm loss for the generator"""
        # First, G(A) should fake the discriminator
        fake_AB = torch.cat((self.real_A, self.fake_B), 1)
        pred_fake, pred_c = self.netD(fake_AB)
        # pred_fake_l = self.netLC(self.fake_B)

        self.loss_G_GAN = self.criterionGAN(pred_fake, True) * 0.6 + \
                          self.criterionGAN(pred_c, True) * 0.4
        # Second, G(A) = B
        # self.loss_G_norm = self.normloss(self.fake_B, self.real_B) * self.opt.lambda_norm
        self.loss_G_norm = self.criterionMSE(self.fake_B, self.real_B)
        self.loss_G_focal = self.focalloss((self.fake_B + 1.) / 2., self.real_label)
        # print(self.loss_M_L1)
        # combine loss and calculate gradients
        self.loss_G = (self.loss_G_norm + self.loss_G_focal) * 15000 + self.loss_G_GAN + \
                       self.criterionMSE(self.fake_length, self.hpwl) * 0.1
        # print(self.loss_G)
        self.loss_G.backward()

        self.GAN_loss += self.loss_G_GAN.sum().item()

    def optimize_parameters(self):
        self.forward()  # compute fake images: G(A)
        # update D, C, LC and LE
        self.set_requires_grad(self.netD, True)  # enable backprop for D
        # self.set_requires_grad(self.netLC, True)
        self.optimizer_D.zero_grad()  # set D's gradients to zero
        # self.optimizer_LC.zero_grad()
        self.backward_D()  # calculate gradients for D
        self.optimizer_D.step()  # update D's weights
        # self.optimizer_LC.step()
        # update G
        self.set_requires_grad(self.netD, False)  # D requires no gradients when optimizing G
        # self.set_requires_grad(self.netLC, False)
        self.optimizer_G.zero_grad()  # set G's gradients to zero
        self.backward_G()  # calculate graidents for G
        self.optimizer_G.step()  # udpate G's weights
