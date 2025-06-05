import numpy as np
import torch
import torch.nn as nn
import os
import math

class Embedder(nn.Module):
    def __init__(self, degree, d_input, d_model):
        super(Embedder, self).__init__()
        self.conv1d = nn.Conv1d(d_input, d_model, 1)
        self.batch_norm = nn.BatchNorm1d(d_model)

    def forward(self, inputs):
        embeddings = self.conv1d(inputs.permute(0, 2, 1))
        embeddings = self.batch_norm(embeddings).permute(0, 2, 1)
        return embeddings

class Pointer(nn.Module):
    def __init__(self, d_query, d_unit):
        super(Pointer, self).__init__()
        self.tanh = nn.Tanh()
        self.w_l = nn.Linear(d_query, d_unit, bias=False)
        self.v = nn.Parameter(torch.FloatTensor(d_unit), requires_grad=True)
        self.v.data.uniform_(-(1. / math.sqrt(d_unit)), 1. / math.sqrt(d_unit))

    def forward(self, refs, query, mask):
        scores = torch.sum(self.v * self.tanh(refs + self.w_l(query).unsqueeze(1)), -1)
        scores = 10. * self.tanh(scores)
        with torch.no_grad():
            scores[mask] = float('-inf')
        return scores

class Glimpse(nn.Module):
    def __init__(self, d_model, d_unit):
        super(Glimpse, self).__init__()
        self.tanh = nn.Tanh()
        self.conv1d = nn.Conv1d(d_model, d_unit, 1)
        self.v = nn.Parameter(torch.FloatTensor(d_unit), requires_grad=True)
        self.v.data.uniform_(-(1. / math.sqrt(d_unit)), 1. / math.sqrt(d_unit))
        self.softmax = nn.Softmax(dim=-1)

    def forward(self, encs):
        encoded = self.conv1d(encs.permute(0, 2, 1)).permute(0, 2, 1)
        scores = torch.sum(self.v * self.tanh(encoded), -1)
        attention = self.softmax(scores)
        glimpse = attention.unsqueeze(-1) * encs
        glimpse = torch.sum(glimpse, 1)
        return glimpse
