import utils
import argparse
from typing import Union

class Config:
    model_class_with_hetero_graph = ["PreRoutGNN"]

    best_count = 0
    best_count_ema = 0
    
    circuits_json_path = None
    model_type = None
    output_dir = None
    seed = None
    test = None
    lr = None
    lr_decay_rate = None
    dropout = None
    gap_update_lr = None
    finetune = None
    save_prediction = None
    pretrained_checkpoint = None
    gap_save_checkpoints = None
    num_epochs = None
    gap_test = None
    gap_save_train_loss = None
    predict_netdelay = None
    predict_celldelay = None
    predict_slew = None
    groundtruth = None
    checkpoints_dir = None
    prediction_dir = None
    root_dir = None
    data_split_json_path = None
    max_gradient_norm = None
    optimizer = None
    comment = None
    num_level_freq_compents = None
    hidden_dim_cellprop = None
    hidden_dim_netprop = None
    hidden_dim_gcn = None
    device = None
    sub_graph_size = None
    num_pin_location_freq_compents = None
    not_check_datasplit = None
    move_to_cuda_in_advance = None
    max_level = None
    normalize_pin_location = None
    use_graph_autoencoder = None
    graph_autoencoder_checkpoint = None
    finetune_graph_autoencoder = None
    graph_autoencoder_hidden_dim = None
    graph_autoencoder_model_type = None
    graph_autoencoder_n_layers = None
    graph_autoencoder_latent_dim = None
    graph_autoencoder_lr = None
    scale_capacitance = None
    homo_GNN_input_node_dim = None
    homo_GNN_input_edge_dim = None
    homo_GNN_output_node_dim = None
    homo_GNN_hidden_dim = None
    homo_GNN_num_layers = None
    ema = None
    ema_decay = None
    ema_gap = None
    ema_start_epoch = None
    scheduler = None


    @staticmethod
    def _save(path:str):
        args = {}
        for k in sorted(dir(Config)):
            if not k.startswith("_"):
                args[k] = getattr(Config, k)
        utils.save_json(args, path)

    @staticmethod
    def _display(n:int=80):
        print("*"*n)
        for k in sorted(dir(Config)):
            if not k.startswith("_"):
                print("{:<30}\t{}".format(k, getattr(Config, k)))
        print("*"*n)

    @staticmethod
    def _load(args:Union[str, dict, argparse.Namespace]):
        if isinstance(args, str):
            args = utils.load_json(args)
        elif isinstance(args, dict):
            pass
        elif isinstance(args, argparse.Namespace):
            args = dict(vars(args))
        else:
            raise NotImplementedError(f"type of args is {type(args)}, not implemented")
        for k,v in args.items():
            setattr(Config, k, v)