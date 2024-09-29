import utils
import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--circuits_json_path', type=str, default="circuits_json/circuits-original.json")
    parser.add_argument('--data_split_json_path', type=str, default="circuits_split_json/data_split-original.json")
    
    parser.add_argument('--root_dir', type=str,default='result-graph_autoencoder', help="Root saving directory")

    parser.add_argument('--model_type', type=str, default="DeepGCNII")
    
    parser.add_argument('--num_level_freq_compents', type=int, default=4, help="Number of frequency components of level encoding. If <= 0, disable")
    parser.add_argument('--max_level', type=int, default=200, help="The max_L in get_positional_encoding, level will divide this argument.")
    parser.add_argument('--num_pin_location_freq_compents', type=int, default=0, help="Number of frequency components of pin location encoding. If <= 0, disable")
    
    parser.add_argument('--hidden_dim', type=int, default=32, help='hidden dim for encoder/decoder')
    parser.add_argument('--latent_dim', type=int, default=4, help='latent dim for encoder/decoder latent space')
    
    parser.add_argument('--lr', type=float, default=5e-4)
    parser.add_argument('--lr_decay_rate', type=float, default=1.0)
    parser.add_argument('--gap_update_lr', type=int, default=100000, help="How many epochs to update learning rate")
    
    parser.add_argument('--dropout', type=float, default=0.0)
    parser.add_argument('--sub_graph_size', type=int, default=50000, help="If graph size > this, partition the graph. If test, it will not partition graph")
    parser.add_argument('--max_gradient_norm', type=float, default=-1)
    parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--num_epochs', type=int, default=10000)
    parser.add_argument('--normalize_pin_location', action='store_true', default=False, help='normalize pin location to [0,1]')
    parser.add_argument('--pretrained_checkpoint', type=str, default=None)
    
    parser.add_argument('--weight_KL_divergence', type=float, default=0.0)
    
    
    parser.add_argument('--n_layers', type=int, default=4, help="Number of layers of homo graph model")
    parser.add_argument('--scale_capacitance', type=float, default=1)
    parser.add_argument('--not_check_datasplit', action='store_true', default=False, help='check if data.json is matched with data_split.json')
    parser.add_argument('--move_to_cuda_in_advance', action='store_true', default=True, help="If True, move all data to cuda together in advance.")
    parser.add_argument('--optimizer', type=str, default="Adam")
    parser.add_argument('--seed', type=int, default=3407)
    parser.add_argument('--groundtruth', action='store_true', default=True, help="Use predecessors' ATSlew to train")
    parser.add_argument('--test', action="store_true", default=False, help="If True, only inference without training")
    parser.add_argument('--finetune', action='store_true', default=False)
    parser.add_argument('--save_inference', action='store_true', default=True, help='Save inference results after training or testing')
    parser.add_argument('--gap_save_checkpoints', type=int, default=200)
    parser.add_argument('--gap_test', type=int, default=20)
    parser.add_argument('--gap_save_train_loss', type=int, default=20)
    parser.add_argument('--comment', type=str, default=None, help="Comment string")
    parser.add_argument('--output_dir', type=str, default=None)
    
    

    args = parser.parse_args()
    return args


def post_process_args(args):    
    if args.test:
        args.finetune = False
        
    if not (args.finetune or args.test):
        args.pretrained_checkpoint = None
        
    if args.test:
        args.sub_graph_size = -1
        args.num_epochs = 1
        args.gap_test = 1
    
    if "cuda" not in args.device:
        args.move_to_cuda_in_advance = False
        
        
    if args.output_dir is None or args.output_dir == "":
        args.output_dir = os.path.join(args.root_dir, args.model_type, utils.get_datetime())
    else:
        args.output_dir = os.path.join(args.root_dir, args.model_type, args.output_dir)
    args.checkpoints_dir = os.path.join(args.output_dir, "checkpoints")
    args.prediction_dir = os.path.join(args.output_dir, "pred")

    return args