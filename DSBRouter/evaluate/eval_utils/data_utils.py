import os
import torchvision.transforms as transforms
import torch

IMG_EXTENSIONS = [
    '.jpg', '.JPG', '.jpeg', '.JPEG',
    '.png', '.PNG', '.ppm', '.PPM', '.bmp', '.BMP',
    '.tif', '.TIF', '.tiff', '.TIFF',
]

def is_image_file(filename):
    return any(filename.endswith(extension) for extension in IMG_EXTENSIONS)

def make_dataset(dirs, max_dataset_size=float("inf")):
    images = []
    for dir in dirs:
        assert os.path.isdir(dir), '%s is not a valid directory' % dir

        for root, _, fnames in sorted(os.walk(dir)):
            for fname in fnames:
                if is_image_file(fname):
                    path = os.path.join(root, fname)
                    images.append(path)
    return images[:min(max_dataset_size, len(images))]
    
def get_transform(grayscale=False, convert=True):
    transform_list = []

    if convert:
        transform_list += [transforms.ToTensor()]
        if grayscale:
            transform_list += [transforms.Normalize((0.5,), (0.5,))]
        else:
            # transform_list += [transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))]
            pass
    return transforms.Compose(transform_list)

def route2vertex(route):
    '''
        route.shape = (batch_size, channel, width, height)
        a pixel is a vertex if 
        1. the left and right pixels has only one non-zero pixel;
        2. the up and down pixels has only one non-zero pixel;
        3. the left, right, up, and down pixels are all non-zero.
    '''
    # change values to 0/1
    route = route / 2. + 0.5
    
    route_up = torch.zeros(route.size())
    route_down = torch.zeros(route.size())
    route_left = torch.zeros(route.size())
    route_right = torch.zeros(route.size())
    route_up[:, :, :-1] = route[:, :, 1:]
    route_down[:, :, 1:] = route[:, :, :-1]
    route_left[:, :-1, :] = route[:, 1:, :]
    route_right[:, 1:, :] = route[:, :-1, :]
    
    route = route * torch.sign((route_up + route_down == 1) + \
        (route_left + route_right == 1) + \
        (route_up + route_down + route_right + route_left == 4))
    
    # change values to -1/1
    route = (route - 0.5) * 2
    return route