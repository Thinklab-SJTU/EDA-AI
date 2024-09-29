from trainer import test
from tqdm import tqdm
from os.path import join
import torch
import random
import torch.nn.functional as F
import pandas as pd
from config import Config
import math
import utils

def save_model(model, path, encoder, encoder_path):
    model.cpu()
    torch.save(model.state_dict(), path)
    model.to(Config.device)
    
    if Config.use_graph_autoencoder:
        encoder.cpu()
        torch.save(encoder.state_dict(), encoder_path)
        encoder.to(Config.device)


# save all things into a single dict
def save_model_packed(path, epoch, model, encoder, optimizer):
    res = {
        "epoch": epoch,
        "model": model.cpu().state_dict(),
        "encoder": encoder.cpu().state_dict() if Config.use_graph_autoencoder else {},
        "optimizer": optimizer.state_dict(),
    }
    torch.save(res, path)
    model.to(Config.device)
    if Config.use_graph_autoencoder:
        encoder.to(Config.device)
    
    

def train(model:torch.nn.Module, data_train:dict, data_test:dict, summary_writer, encoder, model_ema, encoder_ema):
    # csv for total process
    df_train = pd.DataFrame()
    df_test  = pd.DataFrame()
    if Config().ema:
        df_train_ema = pd.DataFrame()
        df_test_ema  = pd.DataFrame()
    
    # csv initial testing
    print("\nInitial testing...")
    if Config.ema:
        df_train_current_epoch_ema, df_test_current_epoch_ema = test(model_ema,  data_train, data_test, summary_writer, -1, encoder_ema)
    df_train_current_epoch, df_test_current_epoch = test(model,  data_train, data_test, summary_writer, -1, encoder)
    print("Initial test finished\n")
        
    df_train = pd.concat([df_train, df_train_current_epoch])
    df_test  = pd.concat([df_test, df_test_current_epoch])
    if Config().ema:
        df_train_ema = pd.concat([df_train_ema, df_train_current_epoch_ema])
        df_test_ema  = pd.concat([df_test_ema, df_test_current_epoch_ema])


    # create optimizer
    optimizer_class = getattr(torch.optim, Config.optimizer)
    if Config.finetune_graph_autoencoder:
        optimizer = optimizer_class([
            {
                "params": model.parameters(),
                "lr": Config.lr,
            },
            {
                "params": encoder.parameters(),
                "lr": Config.graph_autoencoder_lr,             
            },
        ])
    else:
        optimizer = optimizer_class(model.parameters(), lr=Config.lr)

    # load pre-trained optimizer
    if Config.finetune:
        checkpoint = torch.load(Config.pretrained_checkpoint)
        if "optimizer" in checkpoint.keys():
            print("\nLoad optimizer from {}".format(Config.pretrained_checkpoint))
            optimizer.load_state_dict(checkpoint["optimizer"])
            print("\nLoad optimizer from {} finished".format(Config.pretrained_checkpoint))
        
    
    # create lr scheduler
    if Config.scheduler == 'exp':
        scheduler = torch.optim.lr_scheduler.ExponentialLR(optimizer, Config.lr_decay_rate, verbose=False)
    elif Config.scheduler == "cos":
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, Config.num_epochs, Config.lr / 10, verbose=False)
    elif Config.scheduler == "cos_warm_restart":
        scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, 64, 2, Config.lr / 10)
    else:
        raise NotImplementedError(f"lr_scheduler = {Config.scheduler} is not implemented")
    
    
    # judge the best model in testing set
    max_score = -math.inf
    if Config().ema:
        max_score_ema = -math.inf

    for e in tqdm(range(1,Config.num_epochs+1)):
        model.train()

        train_loss_tot_at = 0
        train_loss_tot_slew = 0
        train_loss_tot_netdelay = 0
        train_loss_tot_celldelay = 0

        optimizer.zero_grad()
        
        for k, (g, ts) in random.sample(data_train.items(), len(data_train)):
            # do not modify original graph
            with g.local_scope():
                # move data to cuda
                if not Config.move_to_cuda_in_advance:
                    g, ts = utils.to_device(g, ts, Config.device)
                
                # concat graph embedding
                if Config.use_graph_autoencoder:
                    homo_graph = ts['homo']
                    latent = encoder(homo_graph, homo_graph.ndata['nf'])
                    global_embedding = latent.mean(dim=0).expand_as(latent)
                    g.ndata['nf'] = torch.cat([g.ndata['nf'], latent, global_embedding], dim=1)
                    
                pred_netdelay, pred_celldelay, pred_atslew = model(g, ts, groundtruth=Config.groundtruth)

                loss_at = torch.Tensor([0.0]).to(g.device)
                loss_slew = torch.Tensor([0.0]).to(g.device)
                loss_netdelay = torch.Tensor([0.0]).to(g.device)
                loss_celldelay = torch.Tensor([0.0]).to(g.device)

                loss_at = F.mse_loss(pred_atslew[:,0:4], g.ndata['n_atslew'][:,0:4], reduction="none").mean(dim=1)
                train_loss_tot_at += loss_at.mean().item()

                if Config.predict_slew:
                    loss_slew = F.mse_loss(pred_atslew[:,4:8], g.ndata['n_atslew'][:,4:8], reduction="none").mean(dim=1)
                    train_loss_tot_slew += loss_slew.mean().item()

                if Config.predict_netdelay:
                    loss_netdelay = F.mse_loss(pred_netdelay, g.ndata['n_net_delays_log'], reduction="none").mean(dim=1)
                    train_loss_tot_netdelay += loss_netdelay.mean().item()

                if Config.predict_celldelay:
                    loss_celldelay = F.mse_loss(pred_celldelay, g.edges['cell_out'].data['e_cell_delays'])
                    train_loss_tot_celldelay += loss_celldelay.item()
                 
                if Config.sub_graph_size > 0:
                    loss_at = loss_at * g.ndata['valid']
                    if Config().predict_slew:
                        loss_slew = loss_slew * g.ndata['valid']
                    if Config().predict_netdelay:
                        loss_netdelay = loss_netdelay * g.ndata['valid']
                    

                (loss_at.mean() + loss_slew.mean() + loss_netdelay.mean() + loss_celldelay).backward()
                
                
        if Config.max_gradient_norm > 0:
            gradient_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), Config.max_gradient_norm).item()
        
        # update model after all interation in one epoch
        optimizer.step()
        
        
        if Config.ema and e % Config.ema_gap == 0:
            if e >= Config.ema_start_epoch:
                model_ema.update_parameters(model)
                encoder_ema.update_parameters(encoder)
            else:
                model_ema.module.load_state_dict(model.state_dict())
                encoder_ema.module.load_state_dict(encoder.state_dict())
        

        if e % Config.gap_update_lr == 0: 
            scheduler.step()


        if e % Config.gap_save_train_loss == 0:
            summary_writer.add_scalar("train_loss/AT", train_loss_tot_at, e)
            for model_parameters_index in range(len(optimizer.state_dict()['param_groups'])):
                summary_writer.add_scalar(f"learning_rate/{model_parameters_index}", optimizer.state_dict()['param_groups'][model_parameters_index]['lr'], e)
            if Config.max_gradient_norm > 0: summary_writer.add_scalar("grad_norm", gradient_norm, e)
            if Config.predict_slew: summary_writer.add_scalar("train_loss/Slew", train_loss_tot_slew, e)
            if Config.predict_netdelay: summary_writer.add_scalar("train_loss/netdelay", train_loss_tot_netdelay, e)
            if Config.predict_celldelay: summary_writer.add_scalar("train_loss/celldelay", train_loss_tot_celldelay, e)


        if e % Config.gap_save_checkpoints == 0: 
            save_model_packed(join(Config.checkpoints_dir ,f"{e}.pt"), e, model, encoder, optimizer)
            if Config.ema:
                save_model_packed(join(Config.checkpoints_dir ,f"ema_{e}.pt"), e, model_ema.module, encoder_ema.module, optimizer)
        
        
        if e % Config.gap_test == 0:
            df_train_current_epoch, df_test_current_epoch = test(model, data_train, data_test, summary_writer, e, encoder)
            df_train = pd.concat([df_train, df_train_current_epoch])
            df_test  = pd.concat([df_test, df_test_current_epoch])
            df_train.to_csv(join(Config.output_dir, "train.csv"), index=None)
            df_test.to_csv(join(Config.output_dir, "test.csv"), index=None)
            
            if Config().ema:
                df_train_current_epoch_ema, df_test_current_epoch_ema = test(model_ema, data_train, data_test, summary_writer, e, encoder_ema)
                df_train_ema = pd.concat([df_train_ema, df_train_current_epoch_ema])
                df_test_ema  = pd.concat([df_test_ema, df_test_current_epoch_ema])
                df_train_ema.to_csv(join(Config.output_dir, "train_ema.csv"), index=None)
                df_test_ema.to_csv(join(Config.output_dir, "test_ema.csv"), index=None)
            
            # save best model
            curr_score = df_test.loc[df_test.epoch == e, "r2-AT"].mean()
            if curr_score > max_score:
                max_score = curr_score                
                save_model_packed(join(Config.checkpoints_dir ,f"best_{Config.best_count}.pt"), e, model, encoder, optimizer)
                save_model_packed(join(Config.checkpoints_dir ,f"best.pt"), e, model, encoder, optimizer)
                
                Config.best_count += 1
            
            if Config().ema:
                curr_score_ema = df_test_ema.loc[df_test_ema.epoch == e, "r2-AT"].mean()
                if curr_score_ema > max_score_ema:
                    max_score_ema = curr_score_ema                
                    save_model_packed(join(Config.checkpoints_dir ,f"ema_best_{Config.best_count_ema}.pt"), e, model_ema.module, encoder_ema.module, optimizer)
                    save_model_packed(join(Config.checkpoints_dir ,f"ema_best.pt"), e, model_ema.module, encoder_ema.module, optimizer)
                    
                    Config.best_count_ema += 1       
                                 
    if Config().ema:
        return df_train, df_test, df_train_ema, df_test_ema
    else:
        return df_train, df_test, None, None
        