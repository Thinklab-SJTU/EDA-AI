from . import networks
import torch
import os
import torch.nn as nn
from utils import utils
from abc import ABC


class InferenceModel(nn.Module):

    def __init__(self, opt):
        """Initialize the inference model.

        Parameters:
            opt (Option class)-- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        super().__init__()
        assert (not opt.isTrain)
        self.gpu_ids = opt.gpu_ids
        self.device = torch.device('cuda:{}'.format(self.gpu_ids[0])) if self.gpu_ids else torch.device('cpu')
        self.netG = networks.define_G(3, 1, 64, opt.netG, 'batch', True, 'normal', 0.02, self.gpu_ids)
        self.pretrained_path = opt.path

        # assigns the model to self.netG_[suffix] so that it can be loaded
        # please see <BaseModel.load_networks>

    def initialize(self, epoch):
        load_filename = '%s_net_%s.pth' % (epoch, 'G')
        load_path = os.path.join(self.pretrained_path, load_filename)
        state_dict = torch.load(load_path, map_location=str(self.device))
        net = self.netG.module
        net.load_state_dict(state_dict)
        self.netG.eval()

        print('Loading succeeded.')

    def forward(self, x):
        """Run forward pass."""

        return self.netG(x)
