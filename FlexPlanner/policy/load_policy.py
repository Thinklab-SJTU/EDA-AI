import torch
import utils
from .ppo import PPOPolicy

def load_ppo_policy(ppo_policy:PPOPolicy, checkpoint:str, load_optimizer:bool, device:torch.device) -> int:
    last_epoch = 0
    print(f"[INFO] Load checkpoint from {checkpoint}")
    checkpoint = torch.load(checkpoint, map_location=device)
    utils.load_checkpoint_mismatch(ppo_policy, checkpoint["policy"], allow_mismatch=True)
    print("[INFO] Load PPO successfully")

    if load_optimizer:
        ppo_policy.optim.load_state_dict(checkpoint["optimizer"])
        print("[INFO] Load optimizer successfully")

        last_epoch = checkpoint["epoch"]
        print(f"[INFO] Load epoch {last_epoch} successfully")

    return last_epoch
