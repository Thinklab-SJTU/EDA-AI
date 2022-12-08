import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable


class FocalLoss(nn.Module):
    def __init__(self, gamma=2, alpha=None):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.alpha = alpha

    def forward(self, input, target):
        input = input.view(input.size(0), -1)
        target = target.view(target.size(0), -1)

        eps = 1e-5
        out = -1 * (target * torch.log(0.99 * input + 0.005) * ((target - input) ** self.gamma + eps) +
                    (1 - target) * torch.log(0.995 - 0.99 * input) * (input ** self.gamma + eps))

        return out.mean()


class WCELoss(nn.Module):
    def __init__(self, gamma=50, alpha=None):
        super(WCELoss, self).__init__()
        self.gamma = gamma
        self.alpha = alpha

    def forward(self, input, target):
        input = input.view(input.size(0), -1)
        target = target.view(target.size(0), -1)

        eps = 1e-5
        out = -1 * (target * torch.log(0.99 * input + 0.01) * self.gamma +
                    (1 - target) * torch.log(1. - 0.99 * input))

        return out.mean()


class NormLoss(nn.Module):
    def __init__(self, p=4):
        super(NormLoss, self).__init__()
        self.p = p

    def forward(self, input, target):
        return ((input - target) ** self.p).mean()


class HingeLoss(nn.Module):
    def __init__(self):
        super(HingeLoss, self).__init__()

    def forward(self, pred_real, pred_fake=None):
        if pred_fake is not None:
            loss_real = F.relu(1 - pred_real).mean()
            loss_fake = F.relu(1 + pred_fake).mean()
            return loss_real + loss_fake
        else:
            return -pred_real.mean()


class lengthLoss(nn.Module):
    def __init__(self):
        super(lengthLoss, self).__init__()
        self.relu = nn.ReLU()

    def forward(self, input, target):
        input = input.view(input.size(0), -1)

        input = self.relu(input).sum(1)

        out = torch.abs(input - target)

        return out
