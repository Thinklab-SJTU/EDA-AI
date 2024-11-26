import fp_env
from collections import OrderedDict
import torch
import circuit_dataloader
import utils
import os
from tqdm import tqdm
from arguments import get_args
import math
import tianshou
from tianshou.env import DummyVectorEnv
from tianshou.data import VectorReplayBuffer, Batch
from einops import rearrange
from copy import deepcopy

# from tianshou.policy import PPOPolicy
from policy import PPOPolicy, ClipLossCoef, EntropyLossCoef, load_ppo_policy

# from tianshou.trainer import OnpolicyTrainer
from trainer import OnpolicyTrainer

# from tianshou.data import Collector
from collector import Collector, get_statistics

import model
import numpy as np
from utils import TensorboardWriter
import pandas as pd
from collections import defaultdict
import subprocess
import time


args = get_args()
utils.setup_seed(args.seed)
num_grid_x = args.num_grid_x
num_grid_y = args.num_grid_y
result_dir = args.result_dir
fig_dir = os.path.join(result_dir, "fig")
save_checkpoint_dir = os.path.join(result_dir, "checkpoint")
num_env = args.num_env
num_env_test = args.num_env_test
wiremask_bbo = args.wiremask_bbo
device = torch.device(args.device)
episode_per_collect_per_env = args.episode_per_collect_per_env
save_batch_dir = os.path.join(result_dir, "batch") if args.save_batch else None

utils.mkdir(result_dir)
utils.mkdir(fig_dir)
utils.mkdir(save_checkpoint_dir)
utils.mkdir(save_batch_dir) if save_batch_dir is not None else None
utils.save_json(args.__dict__, os.path.join(result_dir, "args.json"))
utils.save_json({"start": utils.get_datetime()}, os.path.join(result_dir, "time.json"))
hostname = subprocess.check_output("hostname", shell=True).decode().strip()
os.system(f"touch {os.path.join(result_dir, hostname)}")
writer = TensorboardWriter(log_dir=result_dir)

# fp_info
fp_info, df_partner = circuit_dataloader.construct_fp_info_func(args.circuit, args.area_util, num_grid_x, num_grid_y, 
                                                    args.num_alignment, args.alignment_rate, args.alignment_sort, args.num_preplaced_module, args.add_virtual_block, args.num_layer)
episode_len = fp_info.movable_block_num
df_partner.to_csv(os.path.join(result_dir, "partner.csv"), index=False)

# construct shared_encoder
shared_encoder_input_channel = 3 # canvas, wiremask, position mask
if args.input_partner_die:
    shared_encoder_input_channel += (args.num_layer - 1)
if args.input_alignment_mask:
    shared_encoder_input_channel += 1

# one more block one the same layer
if args.input_next_block == 0:
    pass
elif args.input_next_block == 1:
    shared_encoder_input_channel += 2 # wiremask and position mask of next block
else:
    raise NotImplementedError("input_next_block = {} is not implemented".format(args.input_next_block))

shared_encoder_cls = getattr(model, args.shared_encoder_cls)
shared_encoder = model.SharedEncoder(shared_encoder_input_channel, args.graph, args.shared_encoder_final_shape, num_grid_x, shared_encoder_cls, episode_len)

# construct ratio_decider
if args.enable_ratio:
    if args.ratio_share_with_critics:
        ratio_decider_input_dim = shared_encoder.last_channel
        ratio_decider_hidden_dim = 64
    else:
        ratio_decider_input_dim = shared_encoder.input_channels
        ratio_decider_hidden_dim = 4
    ratio_decider = model.RatioDecider(ratio_decider_input_dim, args.ratio_range, args.ratio_area_in_dim, 
                                       ratio_decider_hidden_dim, args.ratio_share_with_critics, args.input_next_block)
else:
    ratio_decider = None

# construct layer_decider
if args.async_place:
    if args.async_place_share_with_critics: # use output from shared encoder
        layer_decider_input_dim = shared_encoder.last_channel
        layer_decider_hidden_dim = 64
    else: # use original stacked mask as input
        layer_decider_input_dim = shared_encoder.input_channels
        layer_decider_hidden_dim = 4

    LayerDeciderClass = model.LayerDecider
    layer_decider = LayerDeciderClass(fp_info.num_layer, layer_decider_input_dim, layer_decider_hidden_dim, 
                                       args.async_place_share_with_critics, args.die_embedding, args.async_place_input_sequence, args.input_layer_sequence)
else:
    layer_decider = None


# local encoder
local_encoder = model.LocalEncoder(shared_encoder.input_channels, 1)

deconv_class = getattr(model, args.deconv_class)
deconv_model = deconv_class(num_grid_x, shared_encoder.final_shape)
actor = model.Actor(
    num_grid_x * num_grid_y, shared_encoder, local_encoder, args.actor_update_shared_encoder, deconv_model=deconv_model, layer_decider=layer_decider,
    wiremask_bbo=wiremask_bbo, ratio_decider=ratio_decider, norm_wiremask=args.norm_wiremask, 
    input_partner_die=args.input_partner_die, input_alignment_mask=args.input_alignment_mask, input_next_block=args.input_next_block, use_ready_layers_mask=False, 
    use_alignment_constraint=args.use_alignment_constraint, set_vision_to_zero=args.set_vision_to_zero, set_canvas_to_zero=args.set_canvas_to_zero, 
)

critic = model.Critic(
    episode_len, shared_encoder, args.norm_wiremask, args.input_partner_die, args.input_alignment_mask, 
    args.input_next_block, fp_info.num_layer, args.input_sequence_critic, args.input_die_critic, args.reduced_dim_critic, args.set_vision_to_zero, args.set_canvas_to_zero, 
)


optimizer = torch.optim.Adam(list(actor.parameters()) + list(critic.parameters()), lr=args.lr)


dist_cls = {
    "pos": torch.distributions.Categorical,
    "ratio": model.VanillaNormal if args.enable_ratio else None,
    "layer": torch.distributions.Categorical if args.async_place else None,
}

clip_loss_coef = ClipLossCoef(args.pos_coef, args.ratio_coef, args.async_place_coef) # args for clip loss
entropy_loss_coef = EntropyLossCoef(async_place_coef=args.async_place_entropy_coef) # args for entropy loss
ppo_policy = PPOPolicy(
    actor, critic, optimizer, dist_cls, max_grad_norm=args.max_grad_norm, 
    ent_coef=args.ent_coef, clip_loss_coef=clip_loss_coef, entropy_loss_coef=entropy_loss_coef, 
    layer_decider_update_gap=args.layer_decider_update_gap, save_batch_dir=save_batch_dir, 
    use_last_step_reward_to_replace_other_steps=args.use_last_step_reward_to_replace_other_steps,
    add_last_step_reward_to_other_steps=args.add_last_step_reward_to_other_steps, error_log_dir=os.path.join(result_dir, "error")
)
ppo_policy.to(device)

with open(os.path.join(result_dir, "model.txt"), "w") as f:
    f.write(str(ppo_policy))

# calculate number of parameters
ppo_params = set(ppo_policy.parameters())
num_ppo_params = sum(p.numel() for p in ppo_params)
utils.save_json({"num_ppo_params": num_ppo_params}, os.path.join(result_dir, "num_ppo_params.json"))    


# load checkpoint first and then collect statistics
last_epoch = 0
if args.load_then_collect and args.checkpoint is not None:
    print("[INFO] load_then_collect = {}".format(args.load_then_collect))
    last_epoch = load_ppo_policy(ppo_policy, args.checkpoint, args.load_optimizer, device)


# print num of parameters of shared_encoder, actor and critic
print(f"[INFO] shared_encoder has {sum(p.numel() for p in shared_encoder.parameters())} parameters")
print(f"[INFO] actor has {sum(p.numel() for p in actor.parameters())} parameters")
print(f"[INFO] critic has {sum(p.numel() for p in critic.parameters())} parameters")


# args for reward function
reward_args = fp_env.RewardArgs(args.reward_func, args.reward_weight_hpwl, args.reward_weight_overlap, args.reward_weight_alignment, args.reward_weight_final_hpwl, 
                                )


# collect statistics for normalization
need_sequence_feature = args.async_place_input_sequence is not None or args.input_sequence_critic
need_alignment_mask = args.input_alignment_mask or args.use_alignment_constraint
single_env = fp_env.PlaceEnv(
    fp_info, args.overlap_ratio, args.along_boundary, reward_args, 
    args.ratio_range, args.async_place, device,
    args.place_order_die_by_die, args.input_next_block,
    args.place_order_sorting_method,
    args.graph, args.input_layer_sequence, need_sequence_feature,
    need_alignment_mask,
)

single_env.reset()
df_place_order = single_env.get_place_order_detailed_information()
df_place_order.to_csv(os.path.join(result_dir, "place_order.csv"), index=False)
buffer_size = num_env * episode_per_collect_per_env * episode_len


# collect statistics
if not wiremask_bbo and args.train:
    if args.statistics is None:
        print("[INFO] Collect statistics for normalization")
        collect_envs = DummyVectorEnv([lambda: deepcopy(single_env) for _ in range(num_env_test)])
        statistics_buffer_size = num_env_test * episode_per_collect_per_env * episode_len
        if args.statistics_method in {1,2}:
            statistics_buffer_size *= 4
        statistics_buffer = VectorReplayBuffer(statistics_buffer_size, num_env_test)
        statistics = get_statistics(collect_envs, statistics_buffer, ppo_policy, statistics_buffer_size // episode_len, args.statistics_method)
        del statistics_buffer, collect_envs
    else:
        print("[INFO] Load statistics from", args.statistics)
        statistics = utils.load_json(args.statistics)
    
    single_env.set_hpwl_norm_coef(statistics["hpwl"])
    utils.save_json(statistics, os.path.join(result_dir, "statistics.json"))

# collect statistics and then load checkpoint
if not args.load_then_collect and args.checkpoint is not None:
    print("[INFO] load_then_collect = {}".format(args.load_then_collect))
    last_epoch = load_ppo_policy(ppo_policy, args.checkpoint, args.load_optimizer, device)

# construct env
train_envs = DummyVectorEnv([lambda: deepcopy(single_env) for _ in range(num_env)])
test_envs = DummyVectorEnv([lambda: deepcopy(single_env) for _ in range(num_env_test)])
train_envs.reset()
test_envs.reset()


# buffer for training
buffer = VectorReplayBuffer(buffer_size, num_env)
buffer.reset()


# collector
train_collector = Collector(ppo_policy, train_envs, buffer)
test_collector = Collector(ppo_policy, test_envs)
train_collector.reset()
test_collector.reset()
test_collector.set_fig_dir(fig_dir, args.save_fig, last_epoch)


# trainer
if not wiremask_bbo and args.train:
    onpolicy_trainer = OnpolicyTrainer(
        ppo_policy,
        train_collector,
        test_collector,
        max_epoch=args.max_epoch,
        step_per_epoch=buffer_size, # collected steps per epoch, collect a buffer
        repeat_per_collect=args.repeat_per_collect, # these collected data will be trained for several times
        episode_per_test=num_env_test, # each env will test one episode
        batch_size=args.batch_size,
        episode_per_collect=buffer_size // episode_len, # collect some episodes for each env
        writer=writer,
        fig_dir=fig_dir,
        save_checkpoint_dir=save_checkpoint_dir,
        save_checkpoint_interval=args.save_checkpoint_interval,
        start_epoch=last_epoch,
        show_progress=False,
        verbose=False,
    )
    result = onpolicy_trainer.run()

print("[INFO] Training finished")

# save tensorboard
writer.save_df()

df = pd.DataFrame()
obs, _ = test_envs.reset()
obs = Batch(obs)
act_record = defaultdict(lambda: defaultdict(list))
ppo_policy.eval()
time_start = time.time()
for iter_idx in tqdm(range(episode_len)):
    # wiremask_bbo = True to use BBO method
    with torch.no_grad():
        probs, _ = ppo_policy.actor(obs)

    act = Batch()
    for key in probs.keys():
        act[key] = ppo_policy.dist_fn[key](probs[key]).sample()
        for env_idx in range(num_env_test):
            single_act = act[key][env_idx].item()

            if key == "ratio":
                next_block_ratio = single_act
                ratio_range = args.ratio_range
                # next_block_ratio is in range [-1,1], use linear transformation to [low, high]
                next_block_ratio = (next_block_ratio + 1) / 2 * (ratio_range[1] - ratio_range[0]) + ratio_range[0]
                next_block_ratio = np.clip(next_block_ratio, ratio_range[0], ratio_range[1])
                single_act = next_block_ratio
                
            act_record[env_idx][key].append(single_act)
    act = tianshou.data.to_numpy(act)
    obs_next, rew, terminated, trunc, info = test_envs.step(act)
    obs_next = Batch(obs_next)

    # # save position mask
    env_idx = 0
    fp_info = test_envs.get_env_attr("fp_info", env_idx)[0]
    if args.save_fig > 0:
        utils.save_intermediate_floorplan(os.path.join(fig_dir, "mask-env={}-place_order={:03d}.png".format(env_idx, iter_idx)), obs["block"][env_idx], 
            obs["canvas"][env_idx],
            obs["position_mask"][env_idx], 
            obs["wiremask"][env_idx], 
            obs["alignment_mask"][env_idx] if "alignment_mask" in obs.keys() else None,
            obs["binary_alignment_mask"][env_idx] if "binary_alignment_mask" in obs.keys() else None,
            fp_info,
        )

    # set next obs
    obs = obs_next

    # add to dataframe
    if terminated.sum() > 0:
        time_end = time.time()
        terminated_indices = np.where(terminated)[0]
        info = Batch(info)
        for terminated_idx in terminated_indices:
            df = pd.concat([df, pd.DataFrame([{
                "env_idx": terminated_idx,
                "hpwl": info["hpwl"][terminated_idx],
                "original_hpwl": info["original_hpwl"][terminated_idx],
                "overlap": info["overlap"][terminated_idx],
                "alignment": info["alignment"][terminated_idx],
                "distance_adjacent_terminal": info["distance_adjacent_terminal"][terminated_idx],
                "max_temp": info["max_temp"][terminated_idx],
                "mean_temp": info["mean_temp"][terminated_idx],
                "runtime": time_end - time_start,
            }])], ignore_index=True)

# save dataframe
df.to_csv(os.path.join(result_dir, "final.csv"), index=False)
utils.save_json(act_record, os.path.join(result_dir, "act_record.json"))
if args.save_fig > 0:
    utils.draw_action_record(act_record, os.path.join(result_dir, "act_record.png"))
time = utils.load_json(os.path.join(result_dir, "time.json"))
time["end"] = utils.get_datetime()
utils.save_json(time, os.path.join(result_dir, "time.json"))

# # save canvas
if args.save_fig > 0 or True:
    for env_idx in range(num_env_test):
        fp_info = test_envs.get_env_attr("fp_info", env_idx)[0]
        utils.save_final_floorplan(os.path.join(fig_dir, "final-canvas-env={:03d}.png".format(env_idx)), fp_info)


lock_path = os.path.join(os.path.dirname(__file__), "lock.tmp")
if os.path.exists(lock_path):
    os.remove(lock_path)
print("[INFO] Done")


with open(os.path.join(result_dir, "DONE"), "w") as f:
    f.write("done")
