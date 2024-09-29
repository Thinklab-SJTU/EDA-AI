# train and test graph auto-encoder
from config import ConfigAutoEncoder
from dataset import process_data_homo
from torch.nn import functional as F
from tqdm import tqdm
from os.path import join
from torch.utils.tensorboard import SummaryWriter
from collections import defaultdict
import model as Model
import pandas as pd
import config
import utils
import torch 
import math

def get_model(in_nf, out_nf, positional_encoding=True, in_ef=12):
    # level encoding, 1+2L, 1 is level, 2L is sin, cos of L frequency
    if ConfigAutoEncoder.num_level_freq_compents > 0 and positional_encoding:
        L = ConfigAutoEncoder.num_level_freq_compents
        in_nf = in_nf + 2*L + 1
    
    # pin location encoding, for each pin, 4 location, each location is encoded with 2*L (sin,cos for L different frequency)
    if ConfigAutoEncoder.num_pin_location_freq_compents > 0 and positional_encoding:
        L = ConfigAutoEncoder.num_pin_location_freq_compents
        in_nf = in_nf + 4*2*L
    
    # in data homo, positinal encoding for ef is not included
    model = getattr(Model, ConfigAutoEncoder.model_type)(
        in_nf,
        in_ef,
        out_nf,
        ConfigAutoEncoder.hidden_dim,
        ConfigAutoEncoder.dropout,
        ConfigAutoEncoder.n_layers,
    )
    return model


def save_checkpoint(encoder, decoder, optimizer, path, epoch, device):
    encoder.cpu()
    decoder.cpu()
    torch.save({
        "epoch": epoch,
        "encoder": encoder.state_dict(),
        "decoder": decoder.state_dict(),
        "optimizer": optimizer.state_dict(),
    }, path)
    encoder.to(device)
    decoder.to(device)


def run(
    encoder,
    decoder,
    mode,
    data,
    df,
    epoch,
    optimizer,
    scheduler,
    summary_writer,
):
    train = True if mode == "train" else False
    
    if train:
        encoder.train()
        decoder.train()
    else:
        encoder.eval()
        decoder.eval()
        
    with torch.set_grad_enabled(train):
        epoch_loss = 0.0
        epoch_loss_reconstruction = 0.0
        epoch_loss_KL = 0.0
        
        if train: 
            optimizer.zero_grad()
            
        for circuit_name, (g,_) in data.items():
            node_feat = g.ndata['nf']
            latent_code = encoder(g, node_feat)
            reconstructed_node_feat = decoder(g, latent_code)
            
            loss_reconstruction = F.mse_loss(reconstructed_node_feat, node_feat[:,0:10], reduction="none").mean(dim=1) * g.ndata['valid']
            loss_reconstruction = loss_reconstruction.mean()
            loss = loss_reconstruction
            
            if ConfigAutoEncoder.weight_KL_divergence > 0:
                latent_var = latent_code.var(dim=0)
                latent_mean = latent_code.mean(dim=0)
                loss_KL = latent_var + latent_mean**2 - 1 - torch.log(latent_var)
                loss_KL = loss_KL.mean()
                loss = loss + ConfigAutoEncoder.weight_KL_divergence * loss_KL
            
            
            epoch_loss_reconstruction += loss_reconstruction.item()
            epoch_loss += loss.item()
            if ConfigAutoEncoder.weight_KL_divergence > 0:
                epoch_loss_KL += loss_KL.item()
            
            
            if train: 
                loss.backward()
            else:
                df = pd.concat([df, pd.DataFrame([{
                    "circuit": circuit_name,
                    "split": mode,
                    "epoch": epoch,
                    "MSE": loss.item(),
                    "loss_reconstruction": loss_reconstruction.item(),
                    "loss_KL": loss_KL.item() if ConfigAutoEncoder.weight_KL_divergence > 0 else None,
                    "loss": loss.item(),
                }])])
                
        if train: 
            if ConfigAutoEncoder.max_gradient_norm > 0:
                gradient_norm_encoder = torch.nn.utils.clip_grad_norm_(encoder.parameters(), ConfigAutoEncoder.max_gradient_norm).item()
                gradient_norm_decoder = torch.nn.utils.clip_grad_norm_(decoder.parameters(), ConfigAutoEncoder.max_gradient_norm).item()
            optimizer.step()
        
        if epoch % ConfigAutoEncoder.gap_update_lr == 0 and train: 
            scheduler.step()
        
        if epoch % ConfigAutoEncoder.gap_save_checkpoints == 0 and train:
            save_checkpoint(encoder, decoder, optimizer, join(ConfigAutoEncoder.checkpoints_dir, f"{epoch}.pt"), epoch, ConfigAutoEncoder.device)
        
        if epoch % ConfigAutoEncoder.gap_save_train_loss == 0: 
            summary_writer.add_scalar(f'{mode}/loss', epoch_loss/len(data), epoch)
            summary_writer.add_scalar(f'{mode}/loss_reconstruction', epoch_loss_reconstruction/len(data), epoch)
            if ConfigAutoEncoder.weight_KL_divergence > 0:
                summary_writer.add_scalar(f'{mode}/loss_KL', epoch_loss_KL/len(data), epoch)
                
            if train and ConfigAutoEncoder.max_gradient_norm > 0:
                summary_writer.add_scalar(f'gradient_norm/encoder', gradient_norm_encoder, epoch)
                summary_writer.add_scalar(f'gradient_norm/decoder', gradient_norm_decoder, epoch)
        
    # use reconstruction to judge the best model
    return df, epoch_loss_reconstruction


def inference(
    encoder,
    decoder,
    g,
):  
    nf = g.ndata['nf']
    latent = encoder(g, nf)
    rec_nf = decoder(g, latent)
    return {
        "nf": nf.detach().cpu(),
        "latent": latent.detach().cpu(),
        "rec_nf": rec_nf.detach().cpu(),
    }
    
                

def main(args):
    utils.setup_seed(args.seed)
    
    ConfigAutoEncoder._load(args)
    ConfigAutoEncoder._display()
    utils.mkdir(ConfigAutoEncoder.output_dir)
    utils.mkdir(ConfigAutoEncoder.checkpoints_dir)
    utils.mkdir(ConfigAutoEncoder.prediction_dir)
    
    ConfigAutoEncoder._save(join(ConfigAutoEncoder.output_dir, "args.json"))
    
    # ===================================== create random encoder =====================================
    # encoder = get_model(10, ConfigAutoEncoder.latent_dim)
    # decoder = get_model(ConfigAutoEncoder.latent_dim, 10, positional_encoding=False)
    # torch.save({
    #     "encoder": encoder.state_dict(),
    #     "decoder": decoder.state_dict(),
    #     "comment": "random model",
    # }, join(ConfigAutoEncoder.checkpoints_dir, "random_model.pt"))
    # exit()
    # ===================================== create random encoder =====================================
    
    data_train, data_test = process_data_homo(
        ConfigAutoEncoder.circuits_json_path, 
        False,
        False, 
        False, 
        ConfigAutoEncoder.data_split_json_path, 
        ConfigAutoEncoder.num_level_freq_compents, 
        ConfigAutoEncoder.sub_graph_size,
        ConfigAutoEncoder.num_pin_location_freq_compents,
        ConfigAutoEncoder.not_check_datasplit,
        ConfigAutoEncoder.move_to_cuda_in_advance,
        ConfigAutoEncoder.device,
        ConfigAutoEncoder.max_level,
        ConfigAutoEncoder.normalize_pin_location,
        False,
        ConfigAutoEncoder.scale_capacitance,
    )
    
    utils.save_json({
        "train": sorted(list(data_train.keys())),
        "test" : sorted(list(data_test.keys())),
    }, join(ConfigAutoEncoder.output_dir, "data_split.json"))
    
    
    encoder = get_model(10, ConfigAutoEncoder.latent_dim).to(ConfigAutoEncoder.device)
    decoder = get_model(ConfigAutoEncoder.latent_dim, 10, positional_encoding=False).to(ConfigAutoEncoder.device)
    if ConfigAutoEncoder.test or ConfigAutoEncoder.finetune:
        print(f"\nLoading checkpoints from {ConfigAutoEncoder.pretrained_checkpoint}...")
        checkpoint = torch.load(ConfigAutoEncoder.pretrained_checkpoint)
        if "encoder" in checkpoint.keys():
            encoder.load_state_dict(checkpoint["encoder"])
        else:
            encoder.load_state_dict(checkpoint)
        if "decoder" in checkpoint.keys():
            decoder.load_state_dict(checkpoint["decoder"])
        else:
            print("\n No checkpoints for decoder, ignored")
            
        if ConfigAutoEncoder.test: 
            print("\nOnly test, inferencing...")
        
    df_train, df_test = pd.DataFrame(), pd.DataFrame()
    
    optimizer = getattr(torch.optim, ConfigAutoEncoder.optimizer)([
        {"params": encoder.parameters()},
        {"params": decoder.parameters()},
    ], lr=ConfigAutoEncoder.lr) if not ConfigAutoEncoder.test else None
    if ConfigAutoEncoder.finetune:
        optimizer.load_state_dict(checkpoint["optimizer"])
        
    scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, ConfigAutoEncoder.lr_decay_rate, verbose=False) if not ConfigAutoEncoder.test else None
    summary_writer = SummaryWriter(ConfigAutoEncoder.output_dir)
    
    best_loss = math.inf
    for epoch in tqdm(range(1,ConfigAutoEncoder.num_epochs+1)):
        if not ConfigAutoEncoder.test:
            run(encoder, decoder, "train", data_train, df_train, epoch, optimizer, scheduler, summary_writer,)
            
        if epoch % ConfigAutoEncoder.gap_test == 0:
            df_train, _         = run(encoder, decoder, "test", data_train, df_train, epoch, None, None, summary_writer,)
            df_test, test_loss  = run(encoder, decoder, "test", data_test,  df_test,  epoch, None, None, summary_writer,)
            
            df_train.to_csv(join(ConfigAutoEncoder.output_dir, "train.csv"), index=None)
            df_test.to_csv(join(ConfigAutoEncoder.output_dir, "test.csv"), index=None)
            pd.concat([df_train, df_test]).to_csv(join(ConfigAutoEncoder.output_dir, "all.csv"), index=None)
            
            # only save checkpoints when enable training
            if test_loss < best_loss:
                best_loss = test_loss
                if not ConfigAutoEncoder.test:
                    save_checkpoint(encoder, decoder, optimizer, join(ConfigAutoEncoder.checkpoints_dir, f"best_{ConfigAutoEncoder.best_count}.pt"), epoch, ConfigAutoEncoder.device)
                    save_checkpoint(encoder, decoder, optimizer, join(ConfigAutoEncoder.checkpoints_dir, f"best.pt"), epoch, ConfigAutoEncoder.device)
                    ConfigAutoEncoder.best_count += 1
    
    
    if ConfigAutoEncoder.save_inference:
        with torch.no_grad():
            if not ConfigAutoEncoder.test:
                checkpoint = torch.load(join(ConfigAutoEncoder.checkpoints_dir, f"best.pt"))
                encoder.load_state_dict(checkpoint["encoder"])
                decoder.load_state_dict(checkpoint["decoder"])
            encoder.eval().to(ConfigAutoEncoder.device)
            decoder.eval().to(ConfigAutoEncoder.device)
            res = defaultdict(dict)
            for k, (g, _) in data_train.items():
                res['train'][k] = inference(encoder, decoder, g)
            for k, (g, _) in data_test.items():
                res['test'][k] = inference(encoder, decoder, g)
            torch.save(res, join(ConfigAutoEncoder.prediction_dir, "node_embedding.pt"))
    
    


if __name__ == "__main__":
    args = config.get_args_autoencoder.get_args()
    args = config.get_args_autoencoder.post_process_args(args)
    main(args)