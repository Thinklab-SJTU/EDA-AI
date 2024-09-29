import torch
import torchmetrics

def calc_r2_flatten(target, pred):
    assert torch.is_tensor(pred)
    assert torch.is_tensor(target)
    return torchmetrics.functional.r2_score(
        pred.detach().flatten(),
        target.detach().flatten()
    )


def calc_r2_torch(target, pred):
    assert torch.is_tensor(pred)
    assert torch.is_tensor(target)
    return torchmetrics.functional.r2_score(
        pred.detach(),
        target.detach()
    )


def calc_mse(target, pred):
    assert torch.is_tensor(pred)
    assert torch.is_tensor(target)
    return torch.nn.functional.mse_loss(pred, target).mean()


def calc_mae(target, pred):
    assert torch.is_tensor(pred)
    assert torch.is_tensor(target)
    return torch.nn.functional.l1_loss(pred, target).mean()


def calc_spearman(target, pred):
    assert torch.is_tensor(pred)
    assert torch.is_tensor(target)
    return torchmetrics.functional.spearman_corrcoef(pred, target).mean()

def calc_pearson(target, pred):
    assert torch.is_tensor(pred)
    assert torch.is_tensor(target)
    return torchmetrics.functional.pearson_corrcoef(pred, target).mean()


def my_r2_fn(target, pred, multioutput="uniform_average"):
    y_truth = target
    y_pred = pred
    ss_res = (y_truth - y_pred).pow(2).sum(dim=0)
    print(f"ss_res is {ss_res}")
    ss_tot = (y_truth - y_truth.mean(dim=0, keepdim=True)).pow(2).sum(dim=0)
    print(f"ss_tot is {ss_tot}")
    res = 1 - ss_res / (ss_tot+1e-8)
    if multioutput == "uniform_average":
        return res.mean()
    elif multioutput == "raw_values":
        return res
    else:
        raise NotImplementedError(multioutput)