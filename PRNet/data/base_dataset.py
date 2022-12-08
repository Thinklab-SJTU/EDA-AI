"""This module implements an abstract base class (ABC) 'BaseDataset' for datasets.

It also includes common transformation functions (e.g., get_transform, __scale_width), which can be later used in subclasses.
"""
import random
import numpy as np
import torch.utils.data as data
from PIL import Image
import torchvision.transforms as transforms
from abc import ABC, abstractmethod


class BaseDataset(data.Dataset, ABC):
    """This class is an abstract base class (ABC) for datasets.

    To create a subclass, you need to implement the following four functions:
    -- <__init__>:                      initialize the class, first call BaseDataset.__init__(self, opt).
    -- <__len__>:                       return the size of dataset.
    -- <__getitem__>:                   get a data point.
    -- <modify_commandline_options>:    (optionally) add dataset-specific options and set default options.
    """

    def __init__(self, opt):
        """Initialize the class; save the options in the class

        Parameters:
            opt (Option class)-- stores all the experiment flags; needs to be a subclass of BaseOptions
        """
        self.opt = opt
        self.root = opt.dataroot

    @staticmethod
    def modify_commandline_options(parser, is_train):
        """Add new dataset-specific options, and rewrite default values for existing options.

        Parameters:
            parser          -- original option parser
            is_train (bool) -- whether training phase or test phase. You can use this flag to add training-specific or test-specific options.

        Returns:
            the modified parser.
        """
        return parser

    @abstractmethod
    def __len__(self):
        """Return the total number of images in the dataset."""
        return 0

    @abstractmethod
    def __getitem__(self, index):
        """Return a data point and its metadata information.

        Parameters:
            index - - a random integer for data indexing

        Returns:
            a dictionary of data with their names. It ususally contains the data itself and its metadata information.
        """
        pass


def crop_pos(img):
    x, y = img.size
    min_x, min_y = x, y
    img_array = np.asarray(img, dtype=np.uint8)
    for i in range(x):
        for j in range(y):
            if img_array[i, j] > 0:
                min_x = min(min_x, i)
                min_y = min(min_y, j)

    return max(0, min_x - 1), max(0, min_y - 1)


def get_params(opt, B):
    # x, y = crop_pos(B)

    # flip = random.random() > 0.5

    # return {'crop_pos': (y, x)}
    return {}


def get_transform(opt, params=None, grayscale=False, method=Image.BICUBIC, convert=True):
    transform_list = []
    # if grayscale:
    #     transform_list.append(transforms.Grayscale(1))

    if 'fill_256' in opt.preprocess:
        transform_list.append(transforms.Lambda(lambda img: __fill256(img, grayscale)))

    if 'fill_512' in opt.preprocess:
        transform_list.append(transforms.Lambda(lambda img: __fill512(img, grayscale)))

    if 'crop' in opt.preprocess:
        transform_list.append(transforms.Lambda(lambda img: __crop(img, params['crop_pos'], opt.crop_size)))

    # if opt.preprocess == 'none':
    #     transform_list.append(transforms.Lambda(lambda img: __make_power_2(img, base=4, method=method)))

    # if not opt.no_flip:
    #     if params is None:
    #         transform_list.append(transforms.RandomHorizontalFlip())
    #     elif params['flip']:
    #         transform_list.append(transforms.Lambda(lambda img: __flip(img, params['flip'])))

    if convert:
        transform_list += [transforms.ToTensor()]
        if grayscale:
            transform_list += [transforms.Normalize((0.5,), (0.5,))]
        else:
            transform_list += [transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))]
    return transforms.Compose(transform_list)


def __fill256(img, grayscale):
    w, h = img.size
    if grayscale:
        _img = Image.new('L', (256, 256), 0)
        _img.paste(img, (0, 0, w, h))
    else:
        _img = Image.new('RGB', (256, 256), (0, 0, 0))
        _img.paste(img, (0, 0, w, h))
    return _img


def __fill512(img, grayscale):
    w, h = img.size
    if grayscale:
        _img = Image.new('L', (512, 512), 0)
        _img.paste(img, (0, 0, w, h))
    else:
        _img = Image.new('RGB', (512, 512), (0, 0, 0))
        _img.paste(img, (0, 0, w, h))
    return _img


def __make_power_2(img, base, method=Image.BICUBIC):
    ow, oh = img.size
    h = int(round(oh / base) * base)
    w = int(round(ow / base) * base)
    if h == oh and w == ow:
        return img

    __print_size_warning(ow, oh, w, h)
    return img.resize((w, h), method)


def __flip(img, flip):
    if flip:
        return img.transpose(Image.FLIP_LEFT_RIGHT)
    return img


def __crop(img, pos, size):
    ow, oh = img.size
    x1, y1 = pos
    tw = th = size
    if (ow > tw or oh > th):
        return img.crop((x1, y1, x1 + tw, y1 + th))
    return img


def __print_size_warning(ow, oh, w, h):
    """Print warning information about image size(only print once)"""
    if not hasattr(__print_size_warning, 'has_printed'):
        print("The image size needs to be a multiple of 4. "
              "The loaded image size was (%d, %d), so it was adjusted to "
              "(%d, %d). This adjustment will be done to all images "
              "whose sizes are not multiples of 4" % (ow, oh, w, h))
        __print_size_warning.has_printed = True
