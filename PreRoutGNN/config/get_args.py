import utils
import argparse
import os
from .config import Config

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--circuits_json_path', type=str, default="circuits_json/circuits-original.json")
    parser.add_argument('--data_split_json_path', type=str, default="circuits_split_json/data_split-original.json")
    
    parser.add_argument('--root_dir', type=str,default='debug', help="Root saving directory")
    
    parser.add_argument('--model_type', type=str, default="PreRoutGNN")
    
    parser.add_argument('--num_level_freq_compents', type=int, default=4, help="Number of frequency components of level encoding. If <= 0, disable")
    parser.add_argument('--max_level', type=int, default=200, help="The max_L in get_positional_encoding")
    parser.add_argument('--num_pin_location_freq_compents', type=int, default=0, help="Number of frequency components of pin location encoding. If <= 0, disable")
    
    # lr and scheduler
    parser.add_argument('--lr', type=float, default=5e-4)
    parser.add_argument('--scheduler', type=str, default='exp')
    parser.add_argument('--lr_decay_rate', type=float, default=1.0, help="lr decay for exp")
    parser.add_argument('--gap_update_lr', type=int, default=100000, help="How many epochs to update learning rate")
    
    
    parser.add_argument('--dropout', type=float, default=0.0)
    parser.add_argument('--max_gradient_norm', type=float, default=1e3)
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--num_epochs', type=int, default=30000)
    parser.add_argument('--normalize_pin_location', action='store_true', default=False, help="Normalize pin location to [0,1]")
    
    
    # for pre-trained graph auto-encoder
    parser.add_argument('--use_graph_autoencoder', action='store_true', default=False, help="use pre-trained graph auto-encoder")
    parser.add_argument('--graph_autoencoder_checkpoint', type=str, default=None)
    parser.add_argument('--finetune_graph_autoencoder', action='store_true', default=True)
    parser.add_argument('--graph_autoencoder_lr', type=float, default=5e-4)
    parser.add_argument('--graph_autoencoder_hidden_dim', type=int, default=32)
    parser.add_argument('--graph_autoencoder_model_type', type=str, default="DeepGCNII")
    parser.add_argument('--graph_autoencoder_n_layers', type=int, default=4)
    parser.add_argument('--graph_autoencoder_latent_dim', type=int, default=4, help="latent space's dimension of auto-encoder")
     
    
    # comparison on HOMO graph model
    parser.add_argument('--homo_GNN_input_node_dim', type=int, default=10)
    parser.add_argument('--homo_GNN_input_edge_dim', type=int, default=12)
    parser.add_argument('--homo_GNN_hidden_dim', type=int, default=32)
    parser.add_argument('--homo_GNN_output_node_dim', type=int, default=4)
    parser.add_argument('--homo_GNN_num_layers', type=int, default=8)
    
    
    # EMA
    parser.add_argument('--ema', action='store_true', default=False)
    parser.add_argument('--ema_decay', type=float, default=0.9)
    parser.add_argument('--ema_gap', type=int, default=1)
    parser.add_argument('--ema_start_epoch', type=int, default=1)
    
    
    # other settings
    parser.add_argument('--pretrained_checkpoint', type=str, default=None, help="use this to test or finetune")
    parser.add_argument('--test', action="store_true", default=False, help="If True, only inference without training")
    parser.add_argument('--sub_graph_size', type=int, default=50000, help="If graph size > this, partition the graph. If test, it will not partition graph. If <= 0, disable")
    parser.add_argument('--move_to_cuda_in_advance', type=int, default=1, help="If True, move all data to cuda together in advance. If > 0, enable; else disable.")
    
    
    # seldom reset these settings
    parser.add_argument('--hidden_dim_netprop', type=int, default=128, help="original is 64")
    parser.add_argument('--hidden_dim_cellprop', type=int, default=128, help="original is 64")
    parser.add_argument('--hidden_dim_gcn', type=int, default=64, help="original is 32")
    parser.add_argument('--scale_capacitance', type=float, default=1.0, help="Capacitance in node feature (nf) will be multiplied by this")
    parser.add_argument('--not_check_datasplit', action='store_true', default=False, help='check if data.json is matched with data_split.json')
    parser.add_argument('--predict_slew', action='store_true', default=True)
    parser.add_argument('--predict_netdelay', action='store_true', default=True)
    parser.add_argument('--predict_celldelay', action='store_true', default=True)
    parser.add_argument('--optimizer', type=str, default="Adam")
    parser.add_argument('--seed', type=int, default=3407)
    parser.add_argument('--groundtruth', action='store_true', default=True, help="Use predecessors' ATSlew to train")
    parser.add_argument('--finetune', action='store_true', default=False)
    parser.add_argument('--save_prediction', action='store_true', default=False)
    parser.add_argument('--gap_save_checkpoints', type=int, default=500)
    parser.add_argument('--gap_test', type=int, default=10)
    parser.add_argument('--gap_save_train_loss', type=int, default=10)
    parser.add_argument('--comment', type=str, default=None, help="Comment string")
    parser.add_argument('--output_dir', type=str, default=None)
    
    args = parser.parse_args()
    return args


def post_process_args(args):    
    if args.move_to_cuda_in_advance > 0:
        args.move_to_cuda_in_advance = True
    else:
        args.move_to_cuda_in_advance = False
        
    if "cuda" not in args.device:
        args.move_to_cuda_in_advance = False
    
    
    
    if args.test and args.finetune:
        raise RuntimeError("cannot test and finetune") 
    
    if not args.finetune and not args.test:
        args.pretrained_checkpoint = None
        
    if args.test:
        args.finetune = False
        args.finetune_graph_autoencoder = False
        args.sub_graph_size = -1
    
    
    
    if args.scheduler != "exp":
        args.lr_decay_rate = None
        
    
        
        
    if args.output_dir is None or args.output_dir == "":
        args.output_dir = os.path.join(args.root_dir, args.model_type, utils.get_datetime())
    else:
        args.output_dir = os.path.join(args.root_dir, args.model_type, args.output_dir)
        
    args.checkpoints_dir = os.path.join(args.output_dir, "checkpoints")
    args.prediction_dir = os.path.join(args.output_dir, "pred")
    
    
    
    if args.num_level_freq_compents <= 0:
        args.max_level = None    
    
    
    if args.model_type not in Config.model_class_with_hetero_graph:
        args.predict_slew = False
        args.predict_netdelay = False
        args.predict_celldelay = False
    else:
        args.homo_GNN_input_node_dim = None
        args.homo_GNN_input_edge_dim = None
        args.homo_GNN_output_node_dim = None
        args.homo_GNN_hidden_dim = None
        args.homo_GNN_num_layers = None
    
    
    
    if not args.use_graph_autoencoder:
        args.graph_autoencoder_checkpoint = None
        args.finetune_graph_autoencoder = None
        args.graph_autoencoder_lr = None
        args.graph_autoencoder_hidden_dim = None
        args.graph_autoencoder_model_type = None
        args.graph_autoencoder_n_layers = None
        args.graph_autoencoder_latent_dim = None
    else:
        if not args.finetune_graph_autoencoder:
            args.lr_graph_autoencoder = None
    
    
    if not args.ema:
        args.ema_decay = None
        args.ema_gap = None
        args.ema_start_epoch = None


    return args