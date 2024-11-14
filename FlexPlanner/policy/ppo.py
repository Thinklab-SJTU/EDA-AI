from typing import Any, Dict, List, Optional, Type

import numpy as np
import torch
from torch import nn

from tianshou.data import Batch, ReplayBuffer, to_torch_as, to_torch
# from tianshou.policy import A2CPolicy
from .a2c import A2CPolicy
from tianshou.utils.net.common import ActorCritic

from tqdm import tqdm
from utils import tqdm_config, DummyTqdm, set_grad, set_grad_none
from collections import defaultdict
import os


class ClipLossCoef:
    def __init__(self, pos_coef:float, ratio_coef:float, async_place_coef:float) -> None:
        self.pos = pos_coef
        self.ratio = ratio_coef
        self.layer = async_place_coef
        # self.layer_next = async_place_coef
    
    def __repr__(self) -> str:
        res = ""
        for key in self.__dict__:
            if key.startswith("_"):
                continue
            res += "\t{:24} = {}\n".format(key+"_coef", self.__dict__[key])
        return res


class EntropyLossCoef:
    def __init__(self, pos_coef:float=1.0, ratio_coef:float=1.0, async_place_coef:float=1.0) -> None:
        self.pos = pos_coef
        self.ratio = ratio_coef
        self.layer = async_place_coef
        # self.layer_next = async_place_coef
    
    def __repr__(self) -> str:
        res = ""
        for key in self.__dict__:
            if key.startswith("_"):
                continue
            res += "\t{:24} = {}\n".format(key+"_coef", self.__dict__[key])
        return res



class PPOPolicy(A2CPolicy):
    r"""Implementation of Proximal Policy Optimization. arXiv:1707.06347.

    :param torch.nn.Module actor: the actor network following the rules in
        :class:`~tianshou.policy.BasePolicy`. (s -> logits)
    :param torch.nn.Module critic: the critic network. (s -> V(s))
    :param torch.optim.Optimizer optim: the optimizer for actor and critic network.
    :param dist_fn: distribution class for computing the action.
    :type dist_fn: Type[torch.distributions.Distribution]
    :param float discount_factor: in [0, 1]. Default to 0.99.
    :param float eps_clip: :math:`\epsilon` in :math:`L_{CLIP}` in the original
        paper. Default to 0.2.
    :param float dual_clip: a parameter c mentioned in arXiv:1912.09729 Equ. 5,
        where c > 1 is a constant indicating the lower bound.
        Default to 5.0 (set None if you do not want to use it).
    :param bool value_clip: a parameter mentioned in arXiv:1811.02553v3 Sec. 4.1.
        Default to True.
    :param bool advantage_normalization: whether to do per mini-batch advantage
        normalization. Default to True.
    :param bool recompute_advantage: whether to recompute advantage every update
        repeat according to https://arxiv.org/pdf/2006.05990.pdf Sec. 3.5.
        Default to False.
    :param float vf_coef: weight for value loss. Default to 0.5.
    :param float ent_coef: weight for entropy loss. Default to 0.01.
    :param float max_grad_norm: clipping gradients in back propagation. Default to
        None.
    :param float gae_lambda: in [0, 1], param for Generalized Advantage Estimation.
        Default to 0.95.
    :param bool reward_normalization: normalize estimated values to have std close
        to 1, also normalize the advantage to Normal(0, 1). Default to False.
    :param int max_batchsize: the maximum size of the batch when computing GAE,
        depends on the size of available memory and the memory cost of the model;
        should be as large as possible within the memory constraint. Default to 256.
    :param bool action_scaling: whether to map actions from range [-1, 1] to range
        [action_spaces.low, action_spaces.high]. Default to True.
    :param str action_bound_method: method to bound action to range [-1, 1], can be
        either "clip" (for simply clipping the action), "tanh" (for applying tanh
        squashing) for now, or empty string for no bounding. Default to "clip".
    :param Optional[gym.Space] action_space: env's action space, mandatory if you want
        to use option "action_scaling" or "action_bound_method". Default to None.
    :param lr_scheduler: a learning rate scheduler that adjusts the learning rate in
        optimizer in each policy.update(). Default to None (no lr_scheduler).
    :param bool deterministic_eval: whether to use deterministic action instead of
        stochastic action sampled by the policy. Default to False.

    .. seealso::

        Please refer to :class:`~tianshou.policy.BasePolicy` for more detailed
        explanation.
    """

    def __init__(
        self,
        actor: torch.nn.Module,
        critic: torch.nn.Module,
        optim: torch.optim.Optimizer,
        dist_fn: dict[str, Type[torch.distributions.Distribution]],
        eps_clip: float = 0.2,
        dual_clip: Optional[float] = None,
        value_clip: bool = False,
        advantage_normalization: bool = True,
        recompute_advantage: bool = False,
        clip_loss_coef: ClipLossCoef = None,
        entropy_loss_coef: EntropyLossCoef = None,
        layer_decider_update_gap:int = 1,
        save_batch_dir:str = None,
        error_log_dir:str = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(actor, critic, optim, dist_fn, **kwargs)
        self._eps_clip = eps_clip
        assert dual_clip is None or dual_clip > 1.0, \
            "Dual-clip PPO parameter should greater than 1.0."
        self._dual_clip = dual_clip
        self._value_clip = value_clip
        self._norm_adv = advantage_normalization
        self._recompute_adv = recompute_advantage
        self._actor_critic: ActorCritic
        self.clip_loss_coef = clip_loss_coef
        assert clip_loss_coef is not None, "clip_loss_coef should not be None"
        print("[INFO] clip_loss_coef:\n", clip_loss_coef)
        self.entropy_loss_coef = entropy_loss_coef
        assert entropy_loss_coef is not None, "entropy_loss_coef should not be None"
        print("[INFO] entropy_loss_coef:\n", entropy_loss_coef)
        self.layer_decider_update_gap = layer_decider_update_gap
        print("[INFO] layer_decider_update_gap = {}".format(layer_decider_update_gap))
        self.save_batch_dir = save_batch_dir
        print("[INFO] save_batch_dir = {}".format(save_batch_dir))

        self.num_error = 0
        self.error_log_dir = error_log_dir
        if self.error_log_dir is not None:
            os.makedirs(self.error_log_dir, exist_ok=True)


    def process_fn(
        self, batch: Batch, buffer: ReplayBuffer, indices: np.ndarray
    ) -> Batch:
        if self._recompute_adv:
            # buffer input `buffer` and `indices` to be used in `learn()`.
            self._buffer, self._indices = buffer, indices
        batch = self._compute_returns(batch, buffer, indices)
        batch.act = to_torch_as(batch.act, batch.v_s)
        with torch.no_grad():
            # batch.logp_old = self(batch).dist.log_prob(batch.act)
            # tmp = self(batch)
            # batch.logp_old = {}
            # for key in tmp.dist.keys():
            #     batch.logp_old[key] = tmp.dist[key].log_prob(batch.act[key].to(device=self.device)).cpu()

            batch.logp_old = batch.info.logp_old
        return batch
    

    def clip_loss(self, dist:torch.distributions.Distribution, act:torch.Tensor, logp_old:torch.Tensor, adv:torch.Tensor) -> torch.Tensor:
        ratio = (dist.log_prob(act) - logp_old).exp().float()
        ratio = ratio.reshape(ratio.size(0), -1).transpose(0, 1)
        surr1 = ratio * adv
        surr2 = ratio.clamp(
            1.0 - self._eps_clip, 1.0 + self._eps_clip
        ) * adv
        if self._dual_clip:
            clip1 = torch.min(surr1, surr2)
            clip2 = torch.max(clip1, self._dual_clip * adv)
            clip_loss = -torch.where(adv < 0, clip2, clip1)
        else:
            clip_loss = -torch.min(surr1, surr2)
        
        # mask = torch.zeros_like(clip_loss)
        # clip_loss = torch.where(clip_loss.abs() < 2, clip_loss, mask)

        return clip_loss.mean()
    

    def learn_async_in_sync_fn(self, probs:torch.Tensor, label:torch.LongTensor, next_block_valid:torch.LongTensor) -> torch.Tensor:
        """
        @probs: (batch_size, num_classes)
        @label: (batch_size,)
        """
        acc = (probs.argmax(dim=1) == label)[next_block_valid == 1].float().mean()
        loss_fn = nn.BCELoss(reduction='none')
        probs_true = torch.zeros_like(probs)
        probs_true[torch.arange(probs.size(0)), label] = 1
        loss = loss_fn(probs, probs_true)
        loss = loss[next_block_valid == 1]
        loss = torch.where(loss <= 10.0, loss, 0)
        loss = loss.mean() / 2 # divide by 2 since each sample is used twice
        return loss, acc
    

    def learn(  # type: ignore
        self, batch: Batch, batch_size: int, repeat: int, **kwargs: Any
    ) -> Dict[str, List[float]]:
        epoch = kwargs['epoch']
        _layer_decider = getattr(self._actor_critic.actor, "layer_decider", None)
        if self.save_batch_dir is not None:
            save_batch_path = os.path.join(self.save_batch_dir, "epoch_{:06d}.pt".format(epoch))
            print("[INFO] save batch to {}".format(save_batch_path))
            torch.save(batch, save_batch_path)
        losses, clip_losses, vf_losses, ent_losses = [], [], [], []

        clip_losses_each = defaultdict(list)
        entropy_losses_each = defaultdict(list)
        num_minibatch = 1 if len(batch) < batch_size else len(batch) // batch_size
        tqdm_class = tqdm
        tqdm_class = DummyTqdm
        progress = tqdm_class(total=repeat * num_minibatch, desc="PPO Learning minibatch, repeat={}, num_minibatch={}".format(repeat, num_minibatch), **tqdm_config)
        for step in range(repeat):
            if self._recompute_adv and step > 0:
                batch = self._compute_returns(batch, self._buffer, self._indices)
            for minibatch in batch.split(batch_size, merge_last=True):
                try:
                    # calculate loss for actor
                    dist = self(minibatch).dist
                    minibatch.adv = minibatch.adv.to(device=self.device)
                    if self._norm_adv:
                        mean, std = minibatch.adv.mean(), minibatch.adv.std()
                        minibatch.adv = (minibatch.adv - mean) / (std + self._eps)  # per-batch norm
                    
                    clip_loss = 0
                    minibatch.act = to_torch(minibatch.act, device=self.device)
                    minibatch.logp_old = to_torch(minibatch.logp_old, device=self.device)
                    for key in dist.keys():
                        clip_loss_curr = self.clip_loss(dist[key], minibatch.act[key], minibatch.logp_old[key], minibatch.adv)
                        clip_loss = clip_loss + clip_loss_curr * getattr(self.clip_loss_coef, key)
                        clip_losses_each[key].append(clip_loss_curr.item())


                    # calculate loss for critic
                    value = self.critic(minibatch.obs).flatten()
                    minibatch.returns = minibatch.returns.to(device=self.device)
                    minibatch.v_s = minibatch.v_s.to(device=self.device)
                    if self._value_clip:
                        v_clip = minibatch.v_s + \
                            (value - minibatch.v_s).clamp(-self._eps_clip, self._eps_clip)
                        vf1 = (minibatch.returns - value).pow(2)
                        vf2 = (minibatch.returns - v_clip).pow(2)
                        vf_loss = torch.max(vf1, vf2).mean()
                    else:
                        vf_loss = (minibatch.returns - value).pow(2).mean()


                    # calculate regularization and overall loss
                    # ent_loss = dist.entropy().mean()
                    ent_loss = 0
                    for key in dist.keys():
                        ent_curr = dist[key].entropy().mean()
                        ent_loss = ent_loss + ent_curr * getattr(self.clip_loss_coef, key) * getattr(self.entropy_loss_coef, key)
                        entropy_losses_each[key].append(ent_curr.item())

                    loss = clip_loss + self._weight_vf * vf_loss \
                        - self._weight_ent * ent_loss
                    self.optim.zero_grad()
                    loss.backward()
                    if self._grad_norm:  # clip large gradient
                        nn.utils.clip_grad_norm_(
                            self._actor_critic.parameters(), max_norm=self._grad_norm
                        )
                    
                    if epoch % self.layer_decider_update_gap == 0:
                        pass
                    else:
                        set_grad_none(_layer_decider)

                    self.optim.step()

                    clip_losses.append(clip_loss.item())
                    vf_losses.append(vf_loss.item())
                    ent_losses.append(ent_loss.item())
                    losses.append(loss.item())


                except ValueError as e:
                    print(e.__repr__())
                    self.num_error += 1
                    if self.error_log_dir is not None:
                        error_log_path = os.path.join(self.error_log_dir, "error-epoch={:06d}-error_idx={:06d}.txt".format(epoch, self.num_error))
                        with open(error_log_path, "w") as f:
                            f.write(e.__repr__())
                finally:
                    progress.update(1)
        

        # total loss
        all_loss = {
            "loss": losses,
            "loss/clip": clip_losses,
            "loss/vf": vf_losses,
            "loss/ent": ent_losses,
        }

        for key in clip_losses_each.keys():
            all_loss["loss/clip_"+key] = clip_losses_each[key]
            all_loss["loss/ent_"+key] = entropy_losses_each[key]
        
        return all_loss
