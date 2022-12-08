import copy
import glob
import os
import time
from collections import deque

import gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

import torch.optim as optim

from a2c_ppo_acktr import algo, utils
from a2c_ppo_acktr.arguments import get_test_args
from a2c_ppo_acktr.storage import RolloutStorage
from evaluation import evaluate
from place_env import place_envs


def main():
    args = get_test_args()

    torch.manual_seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)

    if args.cuda and torch.cuda.is_available() and args.cuda_deterministic:
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True

    log_dir = os.path.expanduser(args.log_dir)
    eval_log_dir = log_dir + "_eval"
    utils.cleanup_log_dir(log_dir)
    utils.cleanup_log_dir(eval_log_dir)

    torch.set_num_threads(1)
    device = torch.device("cuda:0" if args.cuda else "cpu")

    if args.task == 'place':
        envs = place_envs()
        actor_critic = torch.load("./trained_models/placement_300.pt")[0]
        actor_critic.to(device)

    agent = algo.PPO(
            actor_critic,
            args.clip_param,
            args.ppo_epoch,
            args.num_mini_batch,
            args.value_loss_coef,
            args.entropy_coef,
            lr=args.lr,
            eps=args.eps,
            max_grad_norm=args.max_grad_norm)

    rollouts = RolloutStorage(args.num_steps, args.num_processes,
                              envs.obs_space, envs.action_space,
                              actor_critic.recurrent_hidden_state_size)
    obs = envs.reset()
    rollouts.obs[0].copy_(envs.transform(obs))
    rollouts.to(device)
    episode_rewards = deque(maxlen=10)

    features = torch.zeros(710, 2)

    for step in range(args.num_steps):
        # Sample actions
        n = len(envs.results)
        with torch.no_grad():
            value, action, action_log_prob, recurrent_hidden_states = actor_critic.act(
                rollouts.obs[step], rollouts.recurrent_hidden_states[step],
                rollouts.masks[step], features, n)

        # Obser reward and next obs
        obs, done, reward = envs.step(action)
        features[n][0] = action // 32
        features[n][1] = action % 32

        if done:
            obs = envs.reset()
            features = torch.zeros(710, 2)
            print(reward)


if __name__ == "__main__":
    main()
