import torch, os, time, numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation as anime
from diffusers import AutoencoderKL
from torchvision.utils import save_image


class BaseVisualizer():
    def __init__(self, args, device, save_path=None):
        self.args = args
        self.device = device

        self.save_path = save_path
    
    def draw(self):
        raise NotImplementedError


class InferenceResultVisualizer(BaseVisualizer):
    def __init__(self, args, device, save_path=None):
        super().__init__(args, device, save_path)
    
    def draw(self, epoch, iters, xs, xrange=(-1, 1), yrange=(-1, 1), subfix=None):
        save_path = os.path.join(self.save_path,
            f'inference' + (subfix if subfix is not None else ''),
        )
        os.makedirs(save_path, exist_ok=True)
        if self.args.exp2d:
            xs = xs.cpu().numpy()
            plt.figure(figsize=(10, 10))
            plt.scatter(xs[:, 0], xs[:, 1], s=1)
            plt.xlim(xrange[0] - 0.1, xrange[1] + 0.1)
            plt.ylim(yrange[0] - 0.1, yrange[1] + 0.1)
            plt.savefig(os.path.join(save_path, f'ep{epoch}_it{iters}.jpg'))
            plt.close()
        else:
            if xs.shape[1] == 4:
                vae = AutoencoderKL.from_pretrained(
                    "stabilityai/stable-diffusion-xl-base-1.0", subfolder="vae", revision=None, variant=None
                )
                vae = vae.to("cuda")
                samples = xs[:min(16, xs.shape[0])]
                with torch.no_grad():
                    samples = vae.decode(samples.to(self.device) / vae.config.scaling_factor)["sample"].cpu().clamp(-1, 1)

                save_image(samples, (os.path.join(save_path, f'ep{epoch}_it{iters}.png')), nrow=4, normalize=True, value_range=(-1, 1))
                torch.cuda.empty_cache()
            else:
                save_image(xs, (os.path.join(save_path, f'ep{epoch}_it{iters}.png')), nrow=4, normalize=True, value_range=(-1, 1))
                torch.cuda.empty_cache()


class TrajectoryVisualizer(BaseVisualizer):
    def __init__(self, args, device, save_path=None, grid=20):
        super().__init__(args, device, save_path)
        self.grid = grid

    def draw(self, *args, **kwargs):
        if self.args.exp2d:
            self.draw_2d(*args, **kwargs)
        else:
            self.draw_img(*args, **kwargs)
    
    def draw_img(self, epoch, iters, xs, subfix=None):
        save_path = os.path.join(self.save_path,
            f'trajectory' + (subfix if subfix is not None else ''),
            f'ep{epoch}_it{iters}'
        )
        os.makedirs(save_path, exist_ok=True)

        for i in range(xs.shape[1]):
            _xs = xs[:, i]
            if _xs.shape[1] == 4:
                vae = AutoencoderKL.from_pretrained(
                    "stabilityai/stable-diffusion-xl-base-1.0", subfolder="vae", revision=None, variant=None
                )
                vae = vae.to("cuda")
                samples = _xs
                with torch.no_grad():
                    samples = vae.decode(samples.to(self.device) / vae.config.scaling_factor)["sample"].cpu().clamp(-1, 1)

                save_image(samples, (os.path.join(save_path, f'{i}.png')), nrow=10, normalize=True, value_range=(-1, 1))
                torch.cuda.empty_cache()
            else:
                save_image(_xs, (os.path.join(save_path, f'{i}.png')), nrow=10, normalize=True, value_range=(-1, 1))
                torch.cuda.empty_cache()

    def draw_2d(self, epoch, iters, model=None, x_0=None, x_1=None, xs=None, xrange=(-1, 1), yrange=(-1, 1), subfix=None):
        save_path = os.path.join(self.save_path,
            f'trajectory' + (subfix if subfix is not None else ''),
            f'ep{epoch}_it{iters}'
        )
        os.makedirs(save_path, exist_ok=True)

        if xs is None:
            if x_0 is None:
                model.eval()
                x = torch.arange(xrange[0], xrange[1] + 1e-6, (xrange[1] - xrange[0]) / self.grid, dtype=torch.float32, device=self.device)
                y = torch.arange(yrange[0], yrange[1] + 1e-6, (yrange[1] - yrange[0]) / self.grid, dtype=torch.float32, device=self.device)
                x, y = torch.meshgrid(x, x, indexing='xy')
                x_0 = torch.stack([x, y], dim=-1).reshape(-1, 2)
            x_1, xs = model.inference(x_0, return_all=True) if self.args.gpus == 1 else model.module.inference(x_0, return_all=True)
            x_1, xs = x_1.cpu().numpy(), xs.cpu().numpy()
        else:
            x_0, x_1, xs = xs[0].cpu().numpy(), xs[-1].cpu().numpy(), xs.cpu().numpy()

        k = 1000

        plt.figure(figsize=(10, 10))
        plt.scatter(x=np.reshape(xs[:, :k, 0], -1), y=np.reshape(xs[:, :k, 1], -1), s=1, cmap='viridis', vmin=0, vmax=1,
            c=np.reshape(np.repeat(np.expand_dims(np.arange(xs.shape[0]), 1), k, axis=1), -1) / xs.shape[0])
        plt.xlim(xrange[0] - 0.1, xrange[1] + 0.1)
        plt.ylim(yrange[0] - 0.1, yrange[1] + 0.1)
        plt.savefig(os.path.join(save_path, 'trajectory.jpg'))
        plt.close()

        self.draw_animation(xs, save_path, xrange=xrange, yrange=yrange)

    def draw_animation(self, xs, save_path, xrange=(-1, 1), yrange=(-1, 1)):
        clamp = lambda x, a, b: int(min(max(x, a), b))
        st = time.perf_counter()
        num_timesteps, batch_size = xs.shape[0], xs.shape[1]
        steps_per_second = clamp(num_timesteps / 100, 1, 10)
        frames_per_second = clamp(num_timesteps / 10, 1, 10)
        num_seconds = num_timesteps / frames_per_second / steps_per_second + 3
        print('plotting point cloud animation ......', end='', flush=True)
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_xlim(xrange[0] - 0.1, xrange[1] + 0.1)
        ax.set_ylim(yrange[0] - 0.1, yrange[1] + 0.1)
        scatter = ax.scatter([], [], s=1, c=[], cmap='viridis', vmin=0, vmax=1)
        def animate(j):
            j = min((j + 1) * steps_per_second, xs.shape[0])
            cc = np.arange(j) / (num_timesteps - 1)
            cc = np.reshape(np.repeat(np.expand_dims(cc, axis=1), batch_size, axis=1), -1)
            scatter.set_offsets(np.reshape(xs[j - 1], (-1, 2)))
            scatter.set_array(cc[-batch_size:])
            return scatter,
        ani = anime.FuncAnimation(fig, animate, frames=int(num_seconds*frames_per_second), interval=1000/frames_per_second, repeat=False, blit=True)
        try:
            ani.save(os.path.join(save_path, 'trajectory.mp4'), writer=anime.FFMpegWriter(fps=frames_per_second, codec='h264'), dpi=100)
        except:
            ani.save(os.path.join(save_path, 'trajectory.gif'), writer=anime.PillowWriter(fps=frames_per_second), dpi=100)
        plt.close(fig)
        print(f' done! ({time.perf_counter() - st:.3f}s)')