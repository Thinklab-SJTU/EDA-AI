import utils
import argparse
from typing import Union

class ConfigAutoEncoder:
    best_count = 0
    circuits_json_path = None
    model_type = None
    n_layers = None

    output_dir = None
    seed = None
    test = None
    lr = None
    lr_decay_rate = None
    dropout = None
    gap_update_lr = None
    
    finetune = None
    pretrained_checkpoint = None
    gap_save_checkpoints = None
    num_epochs = None
    gap_test = None
    gap_save_train_loss = None
    groundtruth = None
    
    checkpoints_dir = None
    prediction_dir = None
    root_dir = None
    data_split_json_path = None
    max_gradient_norm = None
    optimizer = None
    comment = None
    num_level_freq_compents = None
    hidden_dim = None
    device = None
    sub_graph_size = None
    num_pin_location_freq_compents = None
    not_check_datasplit = None
    move_to_cuda_in_advance = None
    latent_dim = None
    save_inference = None
    max_level = None
    normalize_pin_location = None
    scale_capacitance = None
    weight_KL_divergence = None




    @staticmethod
    def _save(path:str):
        args = {}
        for k in sorted(dir(ConfigAutoEncoder)):
            if not k.startswith("_"):
                args[k] = getattr(ConfigAutoEncoder, k)
        utils.save_json(args, path)

    @staticmethod
    def _display(n:int=80):
        print("*"*n)
        for k in sorted(dir(ConfigAutoEncoder)):
            if not k.startswith("_"):
                print("{:<30}\t{}".format(k, getattr(ConfigAutoEncoder, k)))
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
            setattr(ConfigAutoEncoder, k, v)