import torch 

def my_calc_r2(target, pred):
    y_truth = target
    y_pred = pred
    mse = (y_truth - y_pred).pow(2).mean(dim=0)
    var = (y_truth - y_truth.mean(dim=0, keepdim=True)).pow(2).mean(dim=0)
    r2 = 1 - mse / (var+1e-8)
    r2 = r2.mean()
    return {
        "r2": r2,
        "mse": mse,
        "var": var,
    }

def my_calc_r2_flatten(target, pred):
    y_truth = target.flatten()
    y_pred = pred.flatten()
    mse = (y_truth - y_pred).pow(2).mean(dim=0)
    var = (y_truth - y_truth.mean(dim=0, keepdim=True)).pow(2).mean(dim=0)
    r2 = 1 - mse / (var+1e-8)
    r2 = r2.mean()
    return {
        "r2": r2,
        "mse": mse,
        "var": var,
    }