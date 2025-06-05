import torch, os, tqdm, math
from util.dataset import create_data
from util.noiser import create_noiser
from util.model import create_model
from util.logger import Logger
from util.visualize import InferenceResultVisualizer, TrajectoryVisualizer


class Runner():
    def __init__(self, args):
        self.args = args
        self.rank = self.args.local_rank

        self.dsb = 'dsb' in self.args.method

        self.device = torch.device(f'cuda')

        base_steps_per_epoch = 2 ** 7

        self.prior_set, self.prior_sampler, self.prior_loader = create_data(
            self.args.prior, self.args.gpus, dataset_size=base_steps_per_epoch*self.args.batch_size, 
            batch_size=self.args.batch_size)
        self.data_set, self.data_sampler, self.data_loader = create_data(
            self.args.dataset, self.args.gpus, dataset_size=base_steps_per_epoch*self.args.batch_size, 
            batch_size=self.args.batch_size)
        self.prior_iterator, self.data_iterator = iter(self.prior_loader), iter(self.data_loader)

        self.criterion = torch.nn.MSELoss(reduction='none')
        self.scaler = torch.cuda.amp.GradScaler(enabled=self.args.use_amp)

        if not self.args.exp2d:
            self.val_prior_set = None
            self.val_data_set = None
            if self.args.val_prior is not None and self.args.val_dataset is not None:
                self.val_prior_set, _, _ = create_data(
                    self.args.val_prior, self.args.gpus, dataset_size=base_steps_per_epoch*self.args.batch_size, 
                    batch_size=self.args.batch_size)
                self.val_data_set, _, _ = create_data(
                    self.args.val_dataset, self.args.gpus, dataset_size=base_steps_per_epoch*self.args.batch_size, 
                    batch_size=self.args.batch_size)

        if self.dsb:
            assert self.args.training_timesteps == self.args.inference_timesteps

            self.noiser = create_noiser(self.args.noiser, args, self.device)

            self.cache_size = self.cnt = base_steps_per_epoch * self.args.batch_size * 4

            self.backward_model = create_model(self.args.method, self.args, self.device, self.noiser, rank=self.rank, direction='b')
            self.forward_model = create_model(self.args.method, self.args, self.device, self.noiser, rank=self.rank, direction='f')

            if self.rank == 0:
                print(f'Backward Network #Params: {sum([p.numel() for p in self.backward_model.parameters()])}')
                print(f'Forward Network #Params: {sum([p.numel() for p in self.forward_model.parameters()])}')

            self.backward_model = self.backward_model.to(self.device)
            self.forward_model = self.forward_model.to(self.device)

            self.backward_optimizer = torch.optim.AdamW(
                self.backward_model.parameters(), 
                lr=self.args.lr,
                weight_decay=self.args.weight_decay
            )
            self.forward_optimizer = torch.optim.AdamW(
                self.forward_model.parameters(), 
                lr=self.args.lr,
                weight_decay=self.args.weight_decay
            )

            self.model = {'backward': self.backward_model, 'forward': self.forward_model}
            self.optimizer = {'backward': self.backward_optimizer, 'forward': self.forward_optimizer}
            self.direction = 'backward'
        else:
            self.noiser = create_noiser(self.args.noiser, args, self.device)

            self.model = create_model(self.args.method, self.args, self.device, self.noiser, rank=self.rank)
            if self.rank == 0:
                print(f'#Params: {sum([p.numel() for p in self.model.parameters()])}')

            self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.args.lr)

        self.load_ckpt()

        self.save_path = os.path.join('exp', self.args.exp_name)

        if self.args.global_rank == 0:
            self.evaluators = {
                'Inference': InferenceResultVisualizer(self.args, self.device, save_path=self.save_path),
                'Trajectory': TrajectoryVisualizer(self.args, self.device, save_path=self.save_path),
            }
            self.logger = Logger(os.path.join(self.save_path, 'log'), self.noiser.training_timesteps)

    def _next_batch(self):
        try:
            x_0, x_1 = next(self.prior_iterator), next(self.data_iterator)
        except StopIteration:
            self.prior_iterator, self.data_iterator = iter(self.prior_loader), iter(self.data_loader)
            x_0, x_1 = next(self.prior_iterator), next(self.data_iterator)
        x_0, x_1 = x_0.to(self.device), x_1.to(self.device)
        return x_0, x_1

    def next_batch(self, epoch, dsb=False):
        if dsb:
            if self.cnt + self.args.batch_size > self.cache_size:
                self.x_cache, self.gt_cache, self.t_cache = [], [], []
                self.x_0_cache, self.x_1_cache = [], []
                num_cache = math.ceil(self.cache_size / self.args.batch_size / self.noiser.num_timesteps)
                pbar = tqdm.trange(num_cache, desc=f'Cache epoch {epoch} model {self.direction}') if self.rank == 0 else range(num_cache)
                for _ in pbar:
                    x_0, x_1 = self._next_batch()
                    with torch.no_grad():
                        if self.direction == 'backward' and epoch == 0:
                            _x_cache, _gt_cache, _t_cache = self.noiser.trajectory_dsb(x_0, x_1)
                        else:
                            model = self.model['backward' if self.direction == 'forward' else 'forward']
                            model.eval()
                            _x_cache, _gt_cache, _t_cache = model.inference(x_0 if self.direction == 'forward' else x_1)
                    self.x_cache.append(_x_cache)
                    self.gt_cache.append(_gt_cache)
                    self.t_cache.append(_t_cache)
                self.x_cache = torch.cat(self.x_cache, dim=0).cpu()
                self.gt_cache = torch.cat(self.gt_cache, dim=0).cpu()
                self.t_cache = torch.cat(self.t_cache, dim=0).cpu()
                self.cnt = 0
                self.indexs = torch.randperm(self.x_cache.shape[0])
            index = self.indexs[self.cnt:self.cnt + self.args.batch_size]
            self.cnt += self.args.batch_size
            x = self.x_cache[index]
            gt = self.gt_cache[index]
            t = self.t_cache[index]
            x, gt, t = x.to(self.device), gt.to(self.device), t.to(self.device)
            return x, gt, t
        else:
            x_0, x_1 = self._next_batch()
            t = torch.randint(low=0, high=self.noiser.num_timesteps, size=(x_0.shape[0],), dtype=torch.int64, device=self.device)
            return x_0, x_1, t

    def train(self):
        steps_per_epoch = int(len(self.data_loader) * self.args.repeat_per_epoch)
        print('Steps per epoch:', steps_per_epoch)
        self.cache_size = min(self.cache_size, steps_per_epoch * self.args.batch_size)
        print('Cache size:', self.cache_size)
        print('data_loader:', len(self.data_loader))
        for epoch in range(self.args.epochs):
            if epoch < self.args.skip_epochs:
                print(f'Skipping ep{epoch} and evaluate.')
                self.evaluate(epoch, 0, last=True)
                continue
            self.noiser.train()
            if self.dsb:
                self.cnt = self.cache_size
                self.direction = 'backward' if epoch % 2 == 0 else 'forward'
                model, optimizer = self.model[self.direction], self.optimizer[self.direction]
            else:
                model, optimizer = self.model, self.optimizer
            model.train()
            if self.rank == 0:
                pbar = tqdm.tqdm(total=steps_per_epoch)
            ema_loss, ema_loss_w = None, lambda x: min(0.99, x / 10)
            if self.prior_sampler is not None:
                self.prior_sampler.set_epoch(epoch)
                self.data_sampler.set_epoch(epoch)
            for i in range(steps_per_epoch):
                if self.dsb:
                    x_t, gt, t = self.next_batch(epoch, dsb=True)
                else:
                    x_0, x_1, t = self.next_batch(epoch)
                    x_t = self.noiser(x_0, x_1, t)
                    gt = model.target(x_0, x_1, x_t, t)
                optimizer.zero_grad()
                with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=self.args.use_amp):
                    pred = model(x_t, t)
                    raw_loss = self.criterion(pred, gt).mean(dim=-1)
                    loss = raw_loss.mean()
                ema_loss = loss.item() if ema_loss is None else (ema_loss * ema_loss_w(i) + loss.item() * (1 - ema_loss_w(i)))
                self.scaler.scale(loss).backward()
                self.scaler.step(optimizer)
                self.scaler.update()

                if self.rank == 0 and ((i + 1) % self.args.log_interval == 0 or i == steps_per_epoch - 1):
                    if self.args.global_rank == 0:
                        self.logger.log_step(i, epoch * steps_per_epoch + i, epoch, ema_loss, raw_loss, t, dsb=self.dsb, direction=self.direction if self.dsb else '')
                    desc = f'epoch {epoch}, direction {self.direction}, iteration {i}, loss {ema_loss:.14f}' if self.dsb else f'epoch {epoch}, iteration {i}, loss {ema_loss:.14f}'
                    pbar.set_description(desc, refresh=False)
                    pbar.update(i + 1 - pbar.n)
                if self.args.evaluate_interval >= 0 and (i + 1) % self.args.evaluate_interval == 0 and i != steps_per_epoch - 1:
                    self.evaluate(epoch, i + 1)
                    model.train()

            self.evaluate(epoch, steps_per_epoch, last=True)
            if self.rank == 0:
                pbar.close()

    def evaluate(self, epoch, iters, last=False):
        with torch.no_grad():

            if self.dsb:
                self.backward_model.eval()
                self.forward_model.eval()
            else:
                self.model.eval()

            if self.args.global_rank == 0:

                if last:
                    self.save_ckpt(epoch, iters)

                if not self.args.exp2d:
                    x_prior = torch.stack(
                        [self.val_prior_set[_i] for _i in range(16)], dim=0
                    ).to(self.device)
                    x_data = torch.stack(
                        [self.val_data_set[_i] for _i in range(16)], dim=0
                    ).to(self.device)
                else:
                    x_prior, x_data, _ = self.next_batch(epoch)

                if self.dsb:
                    qs = self.backward_model.inference(x_prior, sample=True)[0]
                    if epoch == 0:
                        ps = self.noiser.trajectory_dsb(x_prior, x_data, sample=True)[0]
                    else:
                        ps = self.forward_model.inference(x_data, sample=True)[0]
                else:
                    qs = self.model.inference(x_prior, return_all=True)[1]
                    ps = self.noiser.trajectory(x_prior, x_data)

                x_prior, x_data, qs, ps = x_prior.to(self.device), x_data.to(self.device), qs.to(self.device), ps.to(self.device)

                os.makedirs(os.path.join(self.save_path, 'trajectory'), exist_ok=True)
                torch.save({
                        'x_prior': x_prior,
                        'x_data': x_data,
                        'qs': qs,
                        'ps': ps,
                    }, os.path.join(self.save_path, 'trajectory', f'ep{epoch}_it{iters}.pth')
                )

                if self.dsb:
                    self.evaluators['Inference'].draw(epoch, iters, qs[-1], subfix=f'_q')
                    self.evaluators['Inference'].draw(epoch, iters, ps[-1], subfix=f'_p')
                else:
                    self.evaluators['Inference'].draw(epoch, iters, qs[-1])

                if last:
                    self.evaluators['Trajectory'].draw(epoch, iters, xs=qs, subfix='_q')
                    self.evaluators['Trajectory'].draw(epoch, iters, xs=ps, subfix='_p')

    def save_ckpt(self, epoch, iters):
        os.makedirs(os.path.join(self.save_path, 'ckpt'), exist_ok=True)
        if self.dsb:
            ckpt = {
                'backward_model': self.backward_model.state_dict(),
                'forward_model': self.forward_model.state_dict(),
                'backward_optimizer': self.backward_optimizer.state_dict(),
                'forward_optimizer': self.forward_optimizer.state_dict(),
            }
        else:
            ckpt = {
                'model': self.model.state_dict(),
                'optimizer': self.optimizer.state_dict(),
            }
        torch.save(ckpt, os.path.join(self.save_path, 'ckpt', f'ep{epoch}_it{iters}.pth'))

    def load_ckpt(self):
        def match_ckpt(ckpt):
            _ckpt = {}
            for k, v in ckpt.items():
                if self.args.gpus > 1 and 'module.' not in k:
                    k = k.replace('network.', 'network.module.')
                elif self.args.gpus == 1 and 'module.' in k:
                    k = k.replace('network.module.', 'network.')
                _ckpt[k] = v
            return _ckpt
        if self.args.ckpt is not None:
            ckpt = torch.load(self.args.ckpt, map_location='cpu')
            if self.dsb:
                self.backward_model.load_state_dict(match_ckpt(ckpt['backward_model']), strict=False)
                self.forward_model.load_state_dict(match_ckpt(ckpt['forward_model']), strict=False)
                if "backward_optimizer" in ckpt:
                    self.backward_optimizer.load_state_dict(ckpt['backward_optimizer'])
                    self.forward_optimizer.load_state_dict(ckpt['forward_optimizer'])
            else:
                self.model.load_state_dict(match_ckpt(ckpt['model']), strict=False)
                if "optimizer" in ckpt:
                    self.optimizer.load_state_dict(ckpt['optimizer'])
        else:
            if self.args.backward_ckpt is not None:
                ckpt = torch.load(self.args.backward_ckpt, map_location='cpu')
                if self.dsb:
                    self.backward_model.load_state_dict(match_ckpt(ckpt['backward_model']), strict=False)
                    if "backward_optimizer" in ckpt:
                        self.backward_optimizer.load_state_dict(ckpt['backward_optimizer'])
                else:
                    self.backward_model.load_state_dict(match_ckpt(ckpt['model']), strict=False)
                    if "optimizer" in ckpt:
                        self.backward_optimizer.load_state_dict(ckpt['optimizer'])
            if self.args.forward_ckpt is not None:
                ckpt = torch.load(self.args.forward_ckpt, map_location='cpu')
                if self.dsb:
                    self.forward_model.load_state_dict(match_ckpt(ckpt['forward_model']), strict=False)
                    if "forward_optimizer" in ckpt:
                        self.forward_optimizer.load_state_dict(ckpt['forward_optimizer'])
                else:
                    self.forward_model.load_state_dict(match_ckpt(ckpt['model']), strict=False)
                    if "optimizer" in ckpt:
                        self.forward_optimizer.load_state_dict(ckpt['optimizer'])

