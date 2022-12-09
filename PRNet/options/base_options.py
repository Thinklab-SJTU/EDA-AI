import argparse
import os
from utils import utils
import torch
import models
import data


class BaseOptions:
    """This class defines options used during both training and test time.

    It also implements several helper functions such as parsing, printing, and saving the options.
    It also gathers additional options defined in <modify_commandline_options> functions in both dataset class and model class.
    """

    def __init__(self):
        """Reset the class; indicates the class hasn't been initailized"""
        self.initialized = False

    def initialize(self, parser):
        """Define the common options that are used in both training and test."""
        # basic parameters
        parser.add_argument('--dataroot', type=str, default='./dataset', help='path to images (should have subfolders trainA, trainB, valA, valB, etc)')
        parser.add_argument('--name', type=str, default='experiment_name', help='name of the experiment. It decides where to store samples and models')
        parser.add_argument('--use_wandb', action='store_true', help='use wandb')
        parser.add_argument('--gpu_ids', type=str, default='0', help='gpu ids: e.g. 0  0,1,2, 0,2. use -1 for CPU')
        parser.add_argument('--checkpoints_dir', type=str, default='./checkpoints', help='models are saved here')
        # model parameters
        parser.add_argument('--model', type=str, default='multi', help='chooses which model to use.')
        parser.add_argument('--input_nc', type=int, default=3, help='# of input image channels: 3 for RGB and 1 for grayscale')
        parser.add_argument('--output_nc', type=int, default=1, help='# of output image channels: 3 for RGB and 1 for grayscale')
        parser.add_argument('--ngf', type=int, default=64, help='# of gen filters in the last conv layer')
        parser.add_argument('--ndf', type=int, default=64, help='# of discrim filters in the first conv layer')
        parser.add_argument('--netD', type=str, default='basic', help='specify discriminator architecture [basic | n_layers | pixel]. The basic model is a 70x70 PatchGAN. n_layers allows you to specify the layers in the discriminator')
        parser.add_argument('--netG', type=str, default='resnet_9blocks', help='specify generator architecture [unet_256 | unet_512 | unet_512_sa]')
        parser.add_argument('--n_layers_D', type=int, default=3, help='only used if netD==n_layers')
        parser.add_argument('--norm', type=str, default='batch', help='instance normalization or batch normalization [instance | batch | none]')
        parser.add_argument('--init_type', type=str, default='normal', help='network initialization [normal | xavier | kaiming | orthogonal]')
        parser.add_argument('--init_gain', type=float, default=0.02, help='scaling factor for normal, xavier and orthogonal.')
        parser.add_argument('--no_dropout', action='store_true', help='no dropout for the generator')
        # dataset parameters
        parser.add_argument('--dataset_mode', type=str, default='aligned', help='chooses how datasets are loaded. [aligned]')
        parser.add_argument('--serial_batches', action='store_true', help='if true, takes images in order to make batches, otherwise takes them randomly')
        parser.add_argument('--crop_size', type=int, default=64, help='then crop to this size')
        parser.add_argument('--max_dataset_size', type=int, default=float("inf"), help='Maximum number of samples allowed per dataset. If the dataset directory contains more than max_dataset_size, only a subset is loaded.')
        parser.add_argument('--preprocess', type=str, default='none', help='scaling and cropping of images at load time [fill | none]')
        parser.add_argument('--no_flip', action='store_true', help='if specified, do not flip the images for data augmentation')
        parser.add_argument('--A_folder', type=str, default='A', help='input A path')
        parser.add_argument('--B_folder', type=str, default='B', help='input B path')
        # additional parameters
        parser.add_argument('--epoch', type=str, default='latest', help='which epoch to load? set to latest to use latest cached model')
        parser.add_argument('--load_iter', type=int, default='0', help='which iteration to load? if load_iter > 0, the code will load models by iter_[load_iter]; otherwise, the code will load models by [epoch]')
        parser.add_argument('--verbose', action='store_true', help='if specified, print more debugging information')
        parser.add_argument('--suffix', default='', type=str, help='customized suffix: opt.name = opt.name + suffix: e.g., {model}_{netG}_size{load_size}')
        parser.add_argument('--path', type=str, default='./pretrained')
        parser.add_argument(
        '--algo', default='a2c', help='algorithm to use: a2c | ppo | acktr')
        parser.add_argument(
            '--gail',
            action='store_true',
            default=False,
            help='do imitation learning with gail')
        parser.add_argument(
            '--gail-experts-dir',
            default='./gail_experts',
            help='directory that contains expert demonstrations for gail')
        parser.add_argument(
            '--gail-batch-size',
            type=int,
            default=128,
            help='gail batch size (default: 128)')
        parser.add_argument(
            '--gail-epoch', type=int, default=5, help='gail epochs (default: 5)')
        parser.add_argument(
            '--lr', type=float, default=7e-4, help='learning rate (default: 7e-4)')
        parser.add_argument(
            '--eps',
            type=float,
            default=1e-5,
            help='RMSprop optimizer epsilon (default: 1e-5)')
        parser.add_argument(
            '--alpha',
            type=float,
            default=0.99,
            help='RMSprop optimizer apha (default: 0.99)')
        parser.add_argument(
            '--gamma',
            type=float,
            default=0.99,
            help='discount factor for rewards (default: 0.99)')
        parser.add_argument(
            '--use-gae',
            action='store_true',
            default=False,
            help='use generalized advantage estimation')
        parser.add_argument(
            '--gae-lambda',
            type=float,
            default=0.95,
            help='gae lambda parameter (default: 0.95)')
        parser.add_argument(
            '--entropy-coef',
            type=float,
            default=0.01,
            help='entropy term coefficient (default: 0.01)')
        parser.add_argument(
            '--value-loss-coef',
            type=float,
            default=0.5,
            help='value loss coefficient (default: 0.5)')
        parser.add_argument(
            '--max-grad-norm',
            type=float,
            default=0.5,
            help='max norm of gradients (default: 0.5)')
        parser.add_argument(
            '--seed', type=int, default=1, help='random seed (default: 1)')
        parser.add_argument(
            '--cuda-deterministic',
            action='store_true',
            default=False,
            help="sets flags for determinism when using CUDA (potentially slow!)")
        parser.add_argument(
            '--num-processes',
            type=int,
            default=16,
            help='how many training CPU processes to use (default: 16)')
        parser.add_argument(
            '--num-steps',
            type=int,
            default=5,
            help='number of forward steps in A2C (default: 5)')
        parser.add_argument(
            '--ppo-epoch',
            type=int,
            default=4,
            help='number of ppo epochs (default: 4)')
        parser.add_argument(
            '--num-mini-batch',
            type=int,
            default=32,
            help='number of batches for ppo (default: 32)')
        parser.add_argument(
            '--clip-param',
            type=float,
            default=0.2,
            help='ppo clip parameter (default: 0.2)')
        parser.add_argument(
            '--log-interval',
            type=int,
            default=10,
            help='log interval, one log per n updates (default: 10)')
        parser.add_argument(
            '--save-interval',
            type=int,
            default=100,
            help='save interval, one save per n updates (default: 100)')
        parser.add_argument(
            '--eval-interval',
            type=int,
            default=None,
            help='eval interval, one eval per n updates (default: None)')
        parser.add_argument(
            '--num-env-steps',
            type=int,
            default=10e6,
            help='number of environment steps to train (default: 10e6)')
        parser.add_argument(
            '--task',
            default='place',
            help='environment to train on (default: place)')
        parser.add_argument(
            '--log-dir',
            default='/tmp/gym/',
            help='directory to save agent logs (default: /tmp/gym)')
        parser.add_argument(
            '--save-dir',
            default='./trained_models/',
            help='directory to save agent logs (default: ./trained_models/)')
        parser.add_argument(
            '--no-cuda',
            action='store_true',
            default=False,
            help='disables CUDA training')
        parser.add_argument(
            '--use-proper-time-limits',
            action='store_true',
            default=False,
            help='compute returns taking into account time limits')
        parser.add_argument(
            '--recurrent-policy',
            action='store_true',
            default=False,
            help='use a recurrent policy')
        parser.add_argument(
            '--use-linear-lr-decay',
            action='store_true',
            default=False,
            help='use a linear schedule on the learning rate')
        parser.add_argument(
            '--log-name',
            type=str,
            default='3-1',
            help='name of save files')
        parser.add_argument(
            '--cell-num',
            type=int,
            default=710,
            help='number of macros')
        parser.add_argument(
            '--net-num',
            type=int,
            default=1454,
            help='number of nets')
        parser.add_argument(
            '--base1',
            type=int,
            default=5800,
            help='baseline one of reward')
        parser.add_argument(
            '--base2',
            type=int,
            default=15000,
            help='baseline two of reward')
        parser.add_argument(
            '--rate1',
            type=float,
            default=0.1,
            help='rate one of reward')
        parser.add_argument(
            '--rate2',
            type=float,
            default=0.04,
            help='rate two of reward')
        self.initialized = True
        return parser

    def gather_options(self):
        """Initialize our parser with basic options(only once).
        Add additional model-specific and dataset-specific options.
        These options are defined in the <modify_commandline_options> function
        in model and dataset classes.
        """
        if not self.initialized:  # check if it has been initialized
            parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            parser = self.initialize(parser)

        # get the basic options
        opt, _ = parser.parse_known_args()

        # modify model-related parser options
        model_name = opt.model
        model_option_setter = models.get_option_setter(model_name)
        parser = model_option_setter(parser, self.isTrain)
        opt, _ = parser.parse_known_args()  # parse again with new defaults

        # modify dataset-related parser options
        dataset_name = opt.dataset_mode
        dataset_option_setter = data.get_option_setter(dataset_name)
        parser = dataset_option_setter(parser, self.isTrain)

        # save and return the parser
        self.parser = parser
        return parser.parse_args()

    def print_options(self, opt):
        """Print and save options

        It will print both current options and default values(if different).
        It will save options into a text file / [checkpoints_dir] / opt.txt
        """
        message = ''
        message += '----------------- Options ---------------\n'
        for k, v in sorted(vars(opt).items()):
            comment = ''
            default = self.parser.get_default(k)
            if v != default:
                comment = '\t[default: %s]' % str(default)
            message += '{:>25}: {:<30}{}\n'.format(str(k), str(v), comment)
        message += '----------------- End -------------------'
        print(message)

        # save to the disk
        expr_dir = os.path.join(opt.checkpoints_dir, opt.name)
        utils.mkdirs(expr_dir)
        file_name = os.path.join(expr_dir, '{}_opt.txt'.format(opt.phase))
        with open(file_name, 'wt') as opt_file:
            opt_file.write(message)
            opt_file.write('\n')

    def parse(self):
        """Parse our options, create checkpoints directory suffix, and set up gpu device."""
        opt = self.gather_options()
        opt.isTrain = self.isTrain   # train or test

        # process opt.suffix
        if opt.suffix:
            suffix = ('_' + opt.suffix.format(**vars(opt))) if opt.suffix != '' else ''
            opt.name = opt.name + suffix

        self.print_options(opt)

        # set gpu ids
        str_ids = opt.gpu_ids.split(',')
        opt.gpu_ids = []
        for str_id in str_ids:
            id = int(str_id)
            if id >= 0:
                opt.gpu_ids.append(id)
        if len(opt.gpu_ids) > 0:
            torch.cuda.set_device(opt.gpu_ids[0])

        self.opt = opt
        return self.opt
