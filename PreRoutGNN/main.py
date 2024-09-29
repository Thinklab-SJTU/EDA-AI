import os
import torch
import utils
import dataset
import pandas as pd
import model as Model

from tqdm import tqdm
from os.path import join
from torch.utils.tensorboard import SummaryWriter
from collections import defaultdict
from trainer import train, test
from config import get_args, post_process_args, Config
from model import ExponentialMovingAverage

# construct timing GNN
def get_model(model_type, predict_slew, dropout):
    out_nf = 4
    if predict_slew:
        out_nf += 4

    # 10=original node embedding, 1=node level, 2*L=level embedding
    in_nf = 10
    if Config.num_level_freq_compents > 0:
        in_nf += (2*Config.num_level_freq_compents + 1)
    
    # 4 location, each location has sin & cos, sin or cos has different L frequency
    if Config.num_pin_location_freq_compents > 0: 
        in_nf += (4*2*Config.num_pin_location_freq_compents)
        
    # global graph embedding & node embedding
    if Config.use_graph_autoencoder:
        in_nf += (2*Config.graph_autoencoder_latent_dim)
    
    # 2 is original, 2*2*L = x,y * 2 * frequency
    ef = 2 
    if Config.num_pin_location_freq_compents > 0:
        ef += (2*2*Config.num_pin_location_freq_compents)
    
    model_class = getattr(Model, model_type)
    model = model_class(in_nf, out_nf, dropout)
    
    return model


# construct timing GNN, for homogeneous graph
def get_homo_model(
    model_type,
    homo_GNN_input_node_dim,
    homo_GNN_input_edge_dim,
    homo_GNN_hidden_dim,
    homo_GNN_output_node_dim,
    homo_GNN_num_layers,
    dropout,
):

    if Config.num_level_freq_compents > 0:
        homo_GNN_input_node_dim += (2*Config.num_level_freq_compents + 1)
    if Config.num_pin_location_freq_compents > 0: 
        homo_GNN_input_node_dim += (4*2*Config.num_pin_location_freq_compents)
    if Config.use_graph_autoencoder:
        homo_GNN_input_node_dim += (2*Config.graph_autoencoder_latent_dim)
    
    model_class = getattr(Model.homo, model_type)
    model = model_class(
        homo_GNN_input_node_dim,
        homo_GNN_input_edge_dim,
        homo_GNN_hidden_dim,
        homo_GNN_output_node_dim,
        homo_GNN_num_layers,
        dropout,
    )
    
    return model



# construct graph auto-encoder
def get_graph_encoder(in_nf, out_nf, in_ef=12):
    # level encoding, 1+2L, 1 is level, 2L is sin, cos of L frequency
    if Config.num_level_freq_compents > 0:
        L = Config.num_level_freq_compents
        in_nf = in_nf + 2*L + 1
    
    # pin location encoding, for each pin, 4 location, each location is encoded with 2*L (sin,cos for L different frequency)
    if Config.num_pin_location_freq_compents > 0:
        L = Config.num_pin_location_freq_compents
        in_nf = in_nf + 4*2*L
    
    # in data homo, positinal encoding for ef is not included
    model = getattr(Model, Config.graph_autoencoder_model_type)(
        in_nf,
        in_ef,
        out_nf,
        Config.graph_autoencoder_hidden_dim,
        Config.dropout,
        Config.graph_autoencoder_n_layers,
    )
    return model



def main(args):
    Config._load(args)
    Config._display()
    utils.setup_seed(Config.seed)


    os.makedirs(Config.output_dir, exist_ok=False)
    os.makedirs(Config.prediction_dir, exist_ok=False)
    os.makedirs(Config.checkpoints_dir, exist_ok=False)
    Config._save(join(Config.output_dir, "args.json"))


    # create model
    if Config.model_type in Config.model_class_with_hetero_graph:
        model = get_model(Config.model_type, Config.predict_slew, Config.dropout)
    else:
        model = get_homo_model(Config.model_type, Config.homo_GNN_input_node_dim, Config.homo_GNN_input_edge_dim, Config.homo_GNN_hidden_dim, Config.homo_GNN_output_node_dim, Config.homo_GNN_num_layers, Config.dropout)
    model.to(Config.device)
    utils.write_model_structure(model, join(Config.output_dir, "model.txt"))
    if Config.test or Config.finetune:
        print(f"\nLoading {Config.model_type} from {Config.pretrained_checkpoint}")
        checkpoint = torch.load(Config.pretrained_checkpoint)
        if "model" in checkpoint.keys():
            model.load_state_dict(checkpoint["model"])
        else:
            model.load_state_dict(checkpoint)
            
    
    
    # pretrained graph autoencoder
    if Config.use_graph_autoencoder:
        print(f"\nLoading pre-trained graph auto-encoder from {Config.graph_autoencoder_checkpoint}")
        encoder = get_graph_encoder(10, Config.graph_autoencoder_latent_dim)
        checkpoint = torch.load(Config.graph_autoencoder_checkpoint)
        if "encoder" in checkpoint.keys():
            encoder.load_state_dict(checkpoint["encoder"])
        else:
            encoder.load_state_dict(checkpoint)
        utils.write_model_structure(encoder, join(Config.output_dir, "encoder.txt"))
    else:
        print(f"\nNot use graph auto-encoder, just nn.Identity()")
        encoder = torch.nn.Identity()
    encoder.train().to(Config.device)
    
    
    # EMA
    if Config.ema:
        print(f"\nEMA for GNN and auto-encoder")
        model_ema = ExponentialMovingAverage(model, Config.ema_decay, Config.device)
        encoder_ema = ExponentialMovingAverage(encoder, Config.ema_decay, Config.device)
    else:
        model_ema = None
        encoder_ema = None
    

    print(f"\nLoading dataset...")
    get_data = dataset.process_data_hetero if Config.model_type in Config.model_class_with_hetero_graph else dataset.process_data_homo
    data_train, data_test = get_data(
        Config.circuits_json_path, 
        Config.predict_slew, 
        Config.predict_netdelay, 
        Config.predict_celldelay, 
        Config.data_split_json_path, 
        Config.num_level_freq_compents, 
        Config.sub_graph_size,
        Config.num_pin_location_freq_compents,
        Config.not_check_datasplit,
        Config.move_to_cuda_in_advance,
        Config.device,
        Config.max_level,
        Config.normalize_pin_location,
        Config.use_graph_autoencoder,
        Config.scale_capacitance,
    )
    


    if Config.test:
        print("\nRunning inference...")
        if Config.ema:
            df_train_ema, df_test_ema = test(model_ema, data_train, data_test, None, -1, encoder_ema)
        df_train, df_test = test(model, data_train, data_test, None, -1, encoder)
            
        
    else:
        summary_writer = SummaryWriter(Config.output_dir)
        data_split = {"train": sorted(list(data_train.keys())), "test": sorted(list(data_test.keys())),}
        utils.save_json(data_split, join(Config.output_dir, "dataset_split.json"))

        df_train, df_test, df_train_ema, df_test_ema = train(model, data_train, data_test, summary_writer, encoder, model_ema, encoder_ema)
        
    # save .csv
    df_train['split'] = "train"
    df_test['split'] = "test"
    
    df_train.to_csv(join(Config.output_dir, "train.csv"), index=None)
    df_test.to_csv(join(Config.output_dir, "test.csv"), index=None)
    
    df_all = pd.concat([df_train, df_test])
    df_all['model'] = Config.model_type
    df_all.to_csv(join(Config.output_dir, "all.csv"), index=None)
    df_all.to_csv(join(Config.output_dir, "{}-{}.csv".format(Config.model_type, Config.seed)), index=None)
    
    
    # save .csv for EMA
    if Config().ema:
        df_train_ema['split'] = "train"
        df_test_ema['split'] = "test"
        
        df_train_ema.to_csv(join(Config.output_dir, "train_ema.csv"), index=None)
        df_test_ema.to_csv(join(Config.output_dir, "test_ema.csv"), index=None)
        
        df_all_ema = pd.concat([df_train_ema, df_test_ema])
        df_all_ema['model'] = Config.model_type
        df_all_ema.to_csv(join(Config.output_dir, "all_ema.csv"), index=None)
    
    print("\nFinished, results are stored in {}".format(Config.output_dir))
    
    
    if Config.test:
        tmp = df_all.loc[
            :,
            [
                "split",
                "circuit",
                "mse-AT",
                "r2-AT",
                "flattenR2-AT",
                "model",
            ]
        ].groupby('split').mean(numeric_only=True)
        print(tmp)



if __name__ == '__main__':
    args = post_process_args(get_args())
    main(args)