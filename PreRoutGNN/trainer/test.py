import torch
import time
import pandas as pd
from os.path import join
from collections import defaultdict
from metric import *
from config import Config
import utils
# from metric.ablation_study import my_calc_r2, my_calc_r2_flatten
import dgl
from model import ExponentialMovingAverage

def test_dict(model, data_dict, epoch, prediction:dict, encoder):
    """this function will modify prediction, which is a dict"""
    
    metric_name_to_metric_func = {
        "mse": calc_mse,
        "r2": calc_r2_torch,
        "flattenR2": calc_r2_flatten,
        # "spearman": calc_spearman,
        # "pearson": calc_pearson,
    }

    predict_slew = Config.predict_slew
    predict_netdelay = Config.predict_netdelay
    predict_celldelay = Config.predict_celldelay
    save_prediction = Config.save_prediction
    prediction_dir = Config.prediction_dir

    

    with torch.no_grad():
        metric_record = defaultdict(dict)
        for circuit_name, (g, ts) in data_dict.items():
            if Config.test:
                print(f"Inference for {circuit_name}")
            
            # for GAT, add self loop
            # g = dgl.add_self_loop(g)
            if not Config.move_to_cuda_in_advance: 
                g, ts = utils.to_device(g, ts, Config.device)
                
            with g.local_scope():
                if "cuda" in Config.device:
                    torch.cuda.synchronize()
                time_s = time.time()
                
                # graph auto encoder
                if Config.use_graph_autoencoder:
                    homo_graph = ts['homo']
                    latent = encoder(homo_graph, homo_graph.ndata['nf'])
                    global_embedding = latent.mean(dim=0).expand_as(latent)
                    g.ndata['nf'] = torch.cat([g.ndata['nf'], latent, global_embedding], dim=1)
                    
                netdelay_pred, celldelay_pred, AT_slew_pred = model(g, ts, groundtruth=False)
                
                if "cuda" in Config.device:
                    torch.cuda.synchronize()
                time_t = time.time()
                time_diff = time_t - time_s


                AT_pred = AT_slew_pred[:, 0:4]
                slew_pred = AT_slew_pred[:, 4:8] if predict_slew else None
                netdelay_pred = netdelay_pred if predict_netdelay else None
                celldelay_pred = celldelay_pred if predict_celldelay else None

                AT_truth = g.ndata['n_atslew'][:, 0:4]
                slew_truth = g.ndata['n_atslew'][:, 4:8] if predict_slew else None
                netdelay_truth = g.ndata['n_net_delays_log'] if predict_netdelay else None
                celldelay_truth = g.edges['cell_out'].data['e_cell_delays'] if predict_celldelay else None
                

                metric_record[circuit_name]["time-time"] = time_diff
                metric_record[circuit_name]["num_nodes-num_nodes"] = g.num_nodes()

                
                if save_prediction:
                    prediction[circuit_name]['at'] = AT_pred.detach().cpu()
                    prediction[circuit_name]['graph'] = g.cpu()
                    if predict_slew: prediction[circuit_name]['slew'] = slew_pred.detach().cpu()
                    if predict_netdelay: prediction[circuit_name]['netdelay'] = netdelay_pred.detach().cpu()
                    if predict_celldelay: prediction[circuit_name]['celldelay'] = celldelay_pred.detach().cpu()
                


                for metric_name, metric_fn in metric_name_to_metric_func.items():
                    r2_AT = metric_fn(AT_truth, AT_pred).item()
                    metric_record[circuit_name][f"{metric_name}-AT"] = r2_AT
                         
                        
                    if predict_slew: 
                        r2_slew = metric_fn(slew_truth, slew_pred).item()
                        metric_record[circuit_name][f"{metric_name}-slew"] = r2_slew

                    if predict_netdelay: 
                        r2_netdelay = metric_fn(netdelay_truth, netdelay_pred).item()
                        metric_record[circuit_name][f"{metric_name}-netdelay"] = r2_netdelay

                    if predict_celldelay: 
                        r2_celldelay = metric_fn(celldelay_truth, celldelay_pred).item()
                        metric_record[circuit_name][f"{metric_name}-celldelay"] = r2_celldelay

        
    return metric_record


def create_dataframe_from_dict(d, epoch):
    df = []
    for circuit in sorted(d.keys()):
        line = {
            "circuit": circuit,
            "epoch": epoch,
        }
        for k,v in d[circuit].items():
            line[k] = v
        df.append(line)
    return pd.DataFrame(df)


def test(model, data_train, data_test, summary_writer, epoch, encoder):
    if Config.test:
        pass
    else:
        data_train = {}
        
    with torch.no_grad():
        predict_slew = Config.predict_slew
        predict_netdelay = Config.predict_netdelay
        predict_celldelay = Config.predict_celldelay

        save_prediction = Config.save_prediction
        prediction_dir = Config.prediction_dir
        
        model.eval()
        if Config.use_graph_autoencoder:
            encoder.eval()

        prediction = defaultdict(dict) if save_prediction else None

        r2_train = test_dict(model, data_train, epoch, prediction, encoder)
        r2_test = test_dict(model, data_test, epoch, prediction, encoder)

        df_train = create_dataframe_from_dict(r2_train, epoch)
        df_test = create_dataframe_from_dict(r2_test, epoch)


        if save_prediction: 
            torch.save(prediction, join(prediction_dir, f"prediction_{epoch}.pt"))
        
        # label = r2-AT, such format
        label_names = list(next(iter(r2_test.values())).keys())
        for label in label_names:
            segs = label.split('-')
            if segs[1] in ["time", "num_nodes"]: continue
            if "slew" in segs[1] and not predict_slew: continue
            if "netdelay" in segs[1] and not predict_netdelay: continue
            if "celldelay" in segs[1] and not predict_celldelay: continue
                            
            for mode in ['train', 'test']:
                r2 = 0
                record = r2_train if mode == "train" else r2_test
                for circuit_name in record.keys():
                    r2 += record[circuit_name][label]
                if len(record)>0 and epoch > -1 and summary_writer is not None:
                    ema = "EMA" if isinstance(model, ExponentialMovingAverage) else "original"
                    summary_writer.add_scalar(f"{mode}_{segs[0]}/{segs[1]}-{ema}", r2/(len(record)), epoch)


        return df_train, df_test