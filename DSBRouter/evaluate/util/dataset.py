import torch, math, numpy as np
from sklearn import datasets
from glob import glob
from PIL import Image
from torchvision import transforms
from diffusers import AutoencoderKL


class BaseDataset(torch.utils.data.Dataset):
    def __init__(self, size=10000, **kwargs):
        self.size = size

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        return self.data[idx]


class CheckerboardDataset(BaseDataset):
    def __init__(self, *args, grid_size=4, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid_size = grid_size
        self.checkboard = torch.tensor([[i, j] for i in range(grid_size) for j in range(grid_size) if (i + j) % 2 == 0])

        grid_pos = torch.randint(low=0, high=self.checkboard.shape[0], size=(self.size,), dtype=torch.int64)
        self.data = torch.rand(size=(self.size, 2), dtype=torch.float32) + self.checkboard[grid_pos].float()
        self.data = self.data / self.grid_size * 2 - 1


class GaussianPointsDataset(BaseDataset):
    def __init__(self, *args, points=16, position=None, sigma=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.points = points
        if position is None:
            gap = 1 / math.sqrt(self.points)
            self.position = [[(i + 0.5) * gap, (j + 0.5) * gap] for i in range(int(math.sqrt(self.points))) for j in range(int(math.sqrt(self.points)))]
        else:
            self.position = position
            assert len(self.position) == self.points
        self.position = torch.tensor(self.position, dtype=torch.float32)
        if sigma is None:
            gap = 1 / math.sqrt(self.points)
            self.sigma = gap / 10
        else:
            self.sigma = sigma

        point_pos = torch.randint(low=0, high=self.points, size=(self.size,), dtype=torch.int64)
        point_pos = self.position[point_pos]
        self.data = torch.randn(size=(self.size, 2), dtype=torch.float32) * self.sigma + point_pos
        self.data = self.data * 2 - 1


class StandardGaussianDataset(BaseDataset):
    def __init__(self, *args, sigma=1, **kwargs):
        super().__init__(*args, **kwargs)
        self.sigma = sigma
        self.data = torch.randn(size=(self.size, 2), dtype=torch.float32) * self.sigma


class LatentStandardGaussianDataset(BaseDataset):
    def __init__(self, *args, sigma=1, **kwargs):
        super().__init__(*args, **kwargs)
        self.sigma = sigma
        self.resolution = kwargs.get('resolution', 32)
        self.data = torch.randn(size=(self.size, 4, 32, 32), dtype=torch.float32) * self.sigma


class PixelStandardGaussianDataset(BaseDataset):
    def __init__(self, *args, sigma=1, **kwargs):
        super().__init__(*args, **kwargs)
        self.sigma = sigma
        self.resolution = kwargs.get('resolution', 32)
        self.data = torch.randn(
            size=(self.size, 3, self.resolution, self.resolution), 
            dtype=torch.float32) * self.sigma

    def __getitem__(self, idx):
        x = torch.randn(
            size=(3, self.resolution, self.resolution), 
            dtype=torch.float32) * self.sigma
        return x


class DSBDataset(torch.utils.data.Dataset):
    def __init__(self, npar, data):
        self.size = npar
        self.data = data

        if data == 'mixture':
            init_sample = torch.randn(npar, 2)
            p = init_sample.shape[0]//2
            init_sample[:p,0] = init_sample[:p,0] * 0.1 - 1/2.
            init_sample[p:,0] = init_sample[p:,0] * 0.1 + 1/2.
            init_sample[:p,1] = init_sample[:p,1] * 0.2
            init_sample[p:,1] = init_sample[p:,1] * 0.2

        elif data == 'scurve':
            X, y  = datasets.make_s_curve(n_samples=npar, noise=0.1, random_state=None)
            init_sample = torch.tensor(X)[:,[0,2]]
            scaling_factor = 0.4
            init_sample = (init_sample - init_sample.mean()) / init_sample.std() * scaling_factor

        elif data == 'swiss':
            X, y  = datasets.make_swiss_roll(n_samples=npar, noise=0.7, random_state=None)
            init_sample = torch.tensor(X)[:,[0,2]]
            scaling_factor = 0.45
            init_sample = (init_sample - init_sample.mean()) / init_sample.std() * scaling_factor        

        elif data == 'moon':
            X, y  = datasets.make_moons(n_samples=npar, noise=0.05, random_state=None)
            scaling_factor = 0.4
            init_sample = torch.tensor(X)
            init_sample = (init_sample - init_sample.mean()) / init_sample.std() * scaling_factor

        elif data == 'circle':
            X, y  = datasets.make_circles(n_samples=npar, noise=0.05, random_state=None, factor=.5)
            init_sample = torch.tensor(X) * 0.8

        elif data == 'checker':
            x1 = np.random.rand(npar) * 4 - 2
            x2_ = np.random.rand(npar) - np.random.randint(0, 2, npar) * 2
            x2 = x2_ + (np.floor(x1) % 2)
            x = np.concatenate([x1[:, None], x2[:, None]], 1) * .5
            init_sample = torch.from_numpy(x)

        elif data == 'pinwheel':
            radial_std = 0.3
            tangential_std = 0.1
            num_classes = 7
            num_per_class = math.ceil(npar / num_classes)
            rate = 0.25
            rads = np.linspace(0, 2 * np.pi, num_classes, endpoint=False)

            features = np.random.randn(num_classes*num_per_class, 2) \
                * np.array([radial_std, tangential_std])
            features[:, 0] += 1.
            labels = np.repeat(np.arange(num_classes), num_per_class)

            angles = rads[labels] + rate * np.exp(features[:, 0])
            rotations = np.stack([np.cos(angles), -np.sin(angles), np.sin(angles), np.cos(angles)])
            rotations = np.reshape(rotations.T, (-1, 2, 2))    
            x = .4 * np.random.permutation(np.einsum("ti,tij->tj", features, rotations))
            init_sample = torch.from_numpy(x)

        elif data == '8gaussians':
            scale = 7.
            centers = [(1, 0), (-1, 0), (0, 1), (0, -1), (1. / np.sqrt(2), 1. / np.sqrt(2)),
                    (1. / np.sqrt(2), -1. / np.sqrt(2)), (-1. / np.sqrt(2),
                                                            1. / np.sqrt(2)), (-1. / np.sqrt(2), -1. / np.sqrt(2))]
            centers = [(scale * x, scale * y) for x, y in centers]

            point = np.random.randn(npar, 2) * 0.5
            idx = np.random.randint(8, size=npar)
            center = np.array([centers[i] for i in idx])
            point += center
            dataset = point
            dataset *= 0.1
            init_sample = torch.from_numpy(dataset)

        if data == '6gaussians':
            scale = 7.
            centers = [(0, 1), (0, -1), (np.sqrt(3) / 2., 1. / 2),
                    (np.sqrt(3) / 2., -1. / 2), (-np.sqrt(3) / 2.,
                                                            1. / 2), (-np.sqrt(3) / 2., -1. / 2)]
            centers = [(scale * x, scale * y) for x, y in centers]

            point = np.random.randn(npar, 2) * 0.5
            idx = np.random.randint(6, size=npar)
            center = np.array([centers[i] for i in idx])
            point += center
            dataset = point
            dataset *= 0.1
            init_sample = torch.from_numpy(dataset)

        self.init_sample = init_sample.float()

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        x = self.init_sample[idx]
        return x


# ------------------------------------------------ #
# afhq dataset                                     #
# ------------------------------------------------ #

class EDADataset(torch.utils.data.Dataset):
    def __init__(self,data_path, size=-1, random_flip=False,
                 resolution=64,center_crop=True, device=None):
        self.data_path = data_path
        self.imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        self.size = size if size > 0 else len(self.imgs)
        self.random_flip = random_flip
        self.device = device

        self.transforms = transforms.Compose(
            [
                transforms.Resize(resolution, interpolation=transforms.InterpolationMode.BILINEAR),
                transforms.CenterCrop(resolution) if center_crop else transforms.RandomCrop(resolution),
                transforms.RandomHorizontalFlip() if random_flip else transforms.Lambda(lambda x: x),
                transforms.ToTensor(),
                transforms.Normalize([0.5], [0.5]),
            ]
        )
        # self.vae = AutoencoderKL.from_pretrained(
        #     "/home/zsh/SDSB/pre_train", revision=None, variant=None
        # )
        # self.vae = self.vae.to(self.device)
    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        img = Image.open(self.imgs[idx % len(self.imgs)]).convert("RGB")
        _data = self.transforms(img)
        # with torch.no_grad():
        #     _data = self.vae.encode(_data.unsqueeze(0)).latent_dist.sample().cpu().squeeze(0) * self.vae.config.scaling_factor
        # print(_data.shape)
        return _data


class NthuDataset(torch.utils.data.Dataset):
    def __init__(self,data_path, size=-1, random_flip=False,
                 resolution=64,center_crop=True, device=None):
        self.data_path = data_path
        self.imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        self.size = size if size > 0 else len(self.imgs)
        self.random_flip = random_flip
        self.device = device

        self.transforms = transforms.Compose(
            [
                transforms.Resize(resolution, interpolation=transforms.InterpolationMode.BILINEAR),
                transforms.CenterCrop(resolution) if center_crop else transforms.RandomCrop(resolution),
                transforms.RandomHorizontalFlip() if random_flip else transforms.Lambda(lambda x: x),
                transforms.ToTensor(),
                transforms.Normalize([0.5], [0.5]),
            ]
        )
        # self.vae = AutoencoderKL.from_pretrained(
        #     "/home/zsh/SDSB/pre_train", revision=None, variant=None
        # )
        # self.vae = self.vae.to(self.device)
    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        img = Image.open(self.imgs[idx % len(self.imgs)]).convert("RGB")
        _data = self.transforms(img)
        # with torch.no_grad():
        #     _data = self.vae.encode(_data.unsqueeze(0)).latent_dist.sample().cpu().squeeze(0) * self.vae.config.scaling_factor
        # print(_data.shape)
        return _data


class ImageAFHQDataset(torch.utils.data.Dataset):
    def __init__(self, data_path, size=-1, random_flip=False,
                 resolution=64, center_crop=True, device=None):
        r"""
        Dataset for AFHQ dataset.

        data_path: str
            Path to the dataset, pre-built as a .png file
        """
        self.data_path = data_path
        self.imgs = sorted(glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg"))
        # self.size = size if size > 0 else len(self.imgs)
        # consistent with the size of 'dog' in afhq
        self.size = 5120
        self.random_flip = random_flip
        self.device = device

        self.train_transforms = transforms.Compose(
            [
                transforms.Resize(resolution, interpolation=transforms.InterpolationMode.BILINEAR),
                transforms.CenterCrop(resolution) if center_crop else transforms.RandomCrop(resolution),
                transforms.RandomHorizontalFlip() if random_flip else transforms.Lambda(lambda x: x),
                transforms.ToTensor(),
                transforms.Normalize([0.5], [0.5]),
            ]
        )

        self.vae = AutoencoderKL.from_pretrained(
            "/inspire/hdd/ws-f4d69b29-e0a5-44e6-bd92-acf4de9990f0/public-project/yanjunchi-24040/zsh/SDSB_afhq/cache/models--stabilityai--stable-diffusion-xl-base-1.0", subfolder="vae", revision=None, variant=None
        )
        self.vae = self.vae.to(self.device)

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        img = Image.open(self.imgs[idx % len(self.imgs)])
        img = img.convert("RGB")
        _data = self.train_transforms(img).to(self.device)
        with torch.no_grad():
            _data = self.vae.encode(_data.unsqueeze(0)).latent_dist.sample().cpu().squeeze(0) * self.vae.config.scaling_factor
        return _data


class ImageCelebADataset(torch.utils.data.Dataset):
    def __init__(self, data_path, size=-1, random_flip=False,
                 resolution=64, center_crop=True,):
        r"""
        Dataset for AFHQ dataset.

        data_path: str
            Path to the dataset, pre-built as a .png file
        """
        self.data_path = data_path
        self.imgs = glob(f"{data_path}/*.png") + glob(f"{data_path}/*.jpg")
        # sort the images
        self.imgs.sort()
        self.size = size if size > 0 else len(self.imgs)
        self.random_flip = random_flip
        
        self.train_transforms = transforms.Compose(
            [
                transforms.Resize(resolution, interpolation=transforms.InterpolationMode.BILINEAR),
                transforms.CenterCrop(resolution) if center_crop else transforms.RandomCrop(resolution),
                transforms.RandomHorizontalFlip() if random_flip else transforms.Lambda(lambda x: x),
                transforms.ToTensor(),
                transforms.Normalize([0.5], [0.5]),
            ]
        )

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        img = Image.open(self.imgs[idx % len(self.imgs)])
        img = img.convert("RGB")
        _data = self.train_transforms(img)
        return _data


def create_data(name, gpus=1, dataset_size=2**24, batch_size=2**16, random_flip=False, device=None):
    name = name.lower()
    if 'checkerboard' in name:
        dataset = CheckerboardDataset(grid_size=int(name.split(':')[-1]), size=dataset_size)
    elif 'gaussians' in name and 'dsb' not in name:
        dataset = GaussianPointsDataset(points=int(name.split(':')[-1]), size=dataset_size)
    elif 'pixel-standard' in name:
        dataset = PixelStandardGaussianDataset(size=dataset_size, resolution=64)
    elif 'latent-standard' in name:
        dataset = LatentStandardGaussianDataset(size=dataset_size)
    elif 'standard' in name:
        dataset = StandardGaussianDataset(size=dataset_size)
    elif 'dsb-scurve' in name:
        dataset = DSBDataset(dataset_size, 'scurve')
    elif 'dsb-swiss' in name:
        dataset = DSBDataset(dataset_size, 'swiss')
    elif 'dsb-moon' in name:
        dataset = DSBDataset(dataset_size, 'moon')
    elif 'dsb-circle' in name:
        dataset = DSBDataset(dataset_size, 'circle')
    elif 'dsb-checker' in name:
        dataset = DSBDataset(dataset_size, 'checker')
    elif 'dsb-pinwheel' in name:
        dataset = DSBDataset(dataset_size, 'pinwheel')
    elif 'dsb-8gaussians' in name:
        dataset = DSBDataset(dataset_size, '8gaussians')
    elif 'dsb-6gaussians' in name:
        dataset = DSBDataset(dataset_size, '6gaussians')
    elif 'afhq-cat' in name:
        dataset = ImageAFHQDataset(
            'dataset/afhq/cat', resolution=int(name.split('-')[-1]),
            size=dataset_size, random_flip=random_flip, device=device)
    elif 'afhq-dog' in name:
        dataset = ImageAFHQDataset(
            'dataset/afhq/dog', resolution=int(name.split('-')[-1]),
            size=dataset_size, random_flip=random_flip, device=device)
    elif 'celeba-64' == name:
        dataset = ImageCelebADataset(
            'dataset/celeba64/',
            size=dataset_size, random_flip=random_flip)   
    elif 'eda_prior' == name:
        dataset = EDADataset(
            'dataset/EDA/prior/',
            size=dataset_size, random_flip=random_flip, device=device)
    elif 'eda_data' == name:
        dataset = EDADataset(
            'dataset/EDA/data/',
            size=dataset_size, random_flip=random_flip, device=device)
    else:
        raise NotImplementedError(f"Dataset {name} not implemented")

    sampler = torch.utils.data.DistributedSampler(dataset) if gpus > 1 else None
    dataloader = torch.utils.data.DataLoader(
        dataset, batch_size=batch_size, 
        shuffle=(sampler is None),
        sampler=sampler,
        num_workers=4)

    return dataset, sampler, dataloader


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    dataset, _, dataloader = create_data('dsb-scurve', dataset_size=2**16)
    xs = next(iter(dataloader))
    print(xs, xs.shape)
    plt.figure(figsize=(10, 10))
    plt.scatter(xs[:, 0], xs[:, 1])
    plt.xlim(-1.1, 1.1)
    plt.ylim(-1.1, 1.1)
    plt.savefig('scurve.jpg')
