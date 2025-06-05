import os, torch, argparse, datetime, random, numpy as np
from util.runner import Runner


def seed_everything(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

def ddp_setup():
    torch.distributed.init_process_group(backend='nccl', timeout=datetime.timedelta(hours=2))
    torch.cuda.set_device(int(os.environ['LOCAL_RANK']))
    print(f'Initialized global rank {os.environ["RANK"]}, local rank {os.environ["LOCAL_RANK"]}')

def main(args):
    args.global_rank = int(os.environ['RANK'])
    args.local_rank = int(os.environ['LOCAL_RANK'])
    args.node = int(args.global_rank // args.gpus)

    ddp_setup()

    seed_everything(42 + args.global_rank)

    runner = Runner(args)

    runner.train()

    return


def create_parser():
    argparser = argparse.ArgumentParser()

    argparser.add_argument('--gpus', type=int, default=8, help='number of gpus per node')
    argparser.add_argument('--node', type=int, default=0, help='number of nodes')
    argparser.add_argument('--local_rank', type=int, default=0, help='local rank')
    argparser.add_argument('--global_rank', type=int, default=0, help='global rank')

    argparser.add_argument('--num_sample', type=int, default=128, help='number of samples')
    
    argparser.add_argument('--exp2d', action='store_true', help='set to true for 2d experiments')
    
    argparser.add_argument('--method', type=str, default='dsb', help='method')
    argparser.add_argument('--simplify', action='store_true', help='whether to use simplified DSB')
    argparser.add_argument('--reparam', type=str, default=None, help='whether to use reparameterized DSB, "term" for TR-DSB, "flow" for FR-DSB')
    argparser.add_argument('--noiser', type=str, default='flow', help='noiser type, "flow" noiser for Flow Matching models, "dsb" noiser for DSB models')
    argparser.add_argument('--gamma_type', type=str, default='constant', help='gamma schedule for DSB')
    argparser.add_argument('--training_timesteps', type=int, default=1000, help='training timesteps')
    argparser.add_argument('--inference_timesteps', type=int, default=100, help='inference timesteps')
    
    argparser.add_argument('--network', type=str, default='mlp', help='network architecture to use')
    argparser.add_argument('--lr', type=float, default=1e-3, help='learning rate')
    argparser.add_argument('--weight_decay', type=float, default=0.0, help='weight decay')
    argparser.add_argument('--batch_size', type=int, default=2**7, help='batch size')
    argparser.add_argument('--epochs', type=int, default=1001, help='number of training epochs')
    argparser.add_argument('--skip_epochs', type=int, default=0, help='number of epochs to skip')
    argparser.add_argument('--repeat_per_epoch', type=float, default=4.0, help='training iteration multiplier per epoch')
    argparser.add_argument('--use_amp', action='store_true', help='whether to use mixed-precision training')
    argparser.add_argument('--log_interval', type=int, default=4, help='interval for printing log')
    argparser.add_argument('--evaluate_interval', type=int, default=8192, help='interval for evaluate and save results')

    argparser.add_argument('--prior', type=str, default='standard', help='prior distribution')
    argparser.add_argument('--dataset', type=str, default='checkerboard:4', help='data distribution')
    argparser.add_argument('--val_prior', type=str, default=None, help='prior distribution for evaluation, only available in image experiments')
    argparser.add_argument('--val_dataset', type=str, default=None, help='data distribution for evaluation, only available in image experiments')

    argparser.add_argument('--exp_name', type=str, default='try', help='name of experiment')
    argparser.add_argument('--ckpt', type=str, default=None, help='checkpoint to load')
    argparser.add_argument('--backward_ckpt', type=str, default=None, help='checkpoint to load for backward model')
    argparser.add_argument('--forward_ckpt', type=str, default=None, help='checkpoint to load for forward model')

    return argparser


if __name__ == '__main__':

    argparser = create_parser()
    args = argparser.parse_args()

    if 'dsb' in args.method:
        assert args.training_timesteps == args.inference_timesteps
    
    main(args)
