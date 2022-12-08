from .base_model import BaseModel
from . import networks
import torch
from utils import utils


class TestModel(BaseModel):
    """ This TestModel can be used to generate testing results.
    This model will automatically set '--dataset_mode test'.
    """

    @staticmethod
    def modify_commandline_options(parser, is_train=True):
        """Add new dataset-specific options, and rewrite default values for existing options.

        Parameters:
            parser          -- original option parser
            is_train (bool) -- whether training phase or test phase. You can use this flag to add training-specific or test-specific options.

        Returns:
            the modified parser.

        The model can only be used during test time. It requires '--dataset_mode single'.
        You need to specify the network using the option '--model_suffix'.
        """
        assert not is_train, 'TestModel cannot be used during training time'
        parser.set_defaults(dataset_mode='test')
        parser.add_argument('--model_suffix', type=str, default='',
                            help='In checkpoints_dir, [epoch]_net_G[model_suffix].pth will be loaded as the generator.')

        return parser

    def __init__(self, opt):
        """Initialize the pix2pix class.

        Parameters:
            opt (Option class)-- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        assert (not opt.isTrain)
        BaseModel.__init__(self, opt)
        # specify the training losses you want to print out. The training/test scripts  will call
        # <BaseModel.get_current_losses>
        self.loss_names = []
        # specify the images you want to save/display. The training/test scripts  will call
        # <BaseModel.get_current_visuals>
        self.visual_names = ['real_A', 'real_B', 'fake_B']
        # specify the models you want to save to the disk. The training/test scripts will call
        # <BaseModel.save_networks> and <BaseModel.load_networks>
        # self.model_names = ['G' + opt.model_suffix, 'D' + opt.model_suffix]  # only generator is needed.
        self.model_names = ['G' + opt.model_suffix]
        self.netG = networks.define_G(opt.input_nc, opt.output_nc, opt.ngf, opt.netG,
                                      opt.norm, not opt.no_dropout, opt.init_type, opt.init_gain, self.gpu_ids)

        # assigns the model to self.netG_[suffix] so that it can be loaded
        # please see <BaseModel.load_networks>
        setattr(self, 'netG' + opt.model_suffix, self.netG)  # store netG in self.

        self.D_accu = 0.
        self.C_accu = 0.

    def set_input(self, input):
        """Unpack input data from the dataloader and perform necessary pre-processing steps.

        Parameters:
            input: a dictionary that contains the data itself and its metadata information.

        We need to use 'single_dataset' dataset mode. It only load images from one domain.
        """
        self.real_A = input['A'].to(self.device)
        self.real_B = input['B'].to(self.device)

    def forward(self):
        """Run forward pass."""
        self.fake_B = self.netG(self.real_A)  # G(real)

        self.realB_length = (self.real_B / 2 + 0.5).view(self.real_B.size(0), -1).sum(1).unsqueeze(-1)
        self.fake_length = torch.round(self.fake_B / 2 + 0.5).view(self.fake_B.size(0), -1).sum(1).unsqueeze(-1)

        self.current_correctness = utils.constraint_label(self.real_A.cpu(),
                                                          self.fake_B.cpu()).unsqueeze(-1).to(self.device)
        self.correctness += self.current_correctness.sum().item()
        self.aligned_real_twl += (self.current_correctness * self.realB_length).sum().item()
        self.aligned_fake_twl += (self.current_correctness * self.fake_length).sum().item()

    def optimize_parameters(self):
        """No optimization for test model."""
        pass
