import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import lightning.pytorch as pl

from torch.distributions.categorical import Categorical
import os
import math

class ScaledDotProductAttention(nn.Module):
    ''' Scaled Dot-Product Attention '''

    def __init__(self, temperature, attn_dropout=0.1):
        super().__init__()
        self.temperature = temperature
        self.dropout = nn.Dropout(attn_dropout)

    def forward(self, q, k, v, mask=None):

        attn = torch.matmul(q / self.temperature, k.transpose(2, 3))

        if mask is not None:
            attn = attn.masked_fill(mask == 0, -1e9)

        attn = self.dropout(F.softmax(attn, dim=-1))
        output = torch.matmul(attn, v)

        return output, attn

class MultiHeadAttention(nn.Module):
    ''' Multi-Head Attention module '''

    def __init__(self, n_head, d_model, d_k, d_v, dropout=0.1):
        super().__init__()

        self.d_model = d_model
        self.n_head = n_head
        self.d_k = d_k
        self.d_v = d_v

        self.w_qs = nn.Linear(d_model, n_head * d_k, bias=False)
        self.w_ks = nn.Linear(d_model, n_head * d_k, bias=False)
        self.w_vs = nn.Linear(d_model, n_head * d_v, bias=False)
        self.fc = nn.Linear(n_head * d_v, d_model, bias=False)

        self.attention = ScaledDotProductAttention(temperature=d_k ** 0.5)

        self.dropout = nn.Dropout(dropout)
        self.batch_norm = nn.BatchNorm1d(d_model)


    def forward(self, q, k, v, mask=None):

        d_k, d_v, n_head = self.d_k, self.d_v, self.n_head
        sz_b, len_q, len_k, len_v = q.size(0), q.size(1), k.size(1), v.size(1)

        residual = q

        # Pass through the pre-attention projection: b x lq x (n*dv)
        # Separate different heads: b x lq x n x dv
        q = self.w_qs(q).view(sz_b, len_q, n_head, d_k)
        k = self.w_ks(k).view(sz_b, len_k, n_head, d_k)
        v = self.w_vs(v).view(sz_b, len_v, n_head, d_v)

        # Transpose for attention dot product: b x n x lq x dv
        q, k, v = q.transpose(1, 2), k.transpose(1, 2), v.transpose(1, 2)

        if mask is not None:
            mask = mask.unsqueeze(1)   # For head axis broadcasting.

        q, attn = self.attention(q, k, v, mask=mask)

        # Transpose to move the head dimension back: b x lq x n x dv
        # Combine the last two dimensions to concatenate all the heads together: b x lq x (n*dv)
        q = q.transpose(1, 2).contiguous().view(sz_b, len_q, -1)
        q = self.dropout(self.fc(q))
        q += residual

        seq_len = q.size(1)
        q = self.batch_norm(q.view(-1, self.d_model)).view(-1, seq_len, self.d_model)

        return q, attn

class PositionwiseFeedForward(nn.Module):
    ''' A two-feed-forward-layer module '''

    def __init__(self, d_in, d_hid, dropout=0.1):
        super().__init__()
        self.d_in = d_in
        self.w_1 = nn.Linear(d_in, d_hid) # position-wise
        self.w_2 = nn.Linear(d_hid, d_in) # position-wise
        self.batch_norm = nn.BatchNorm1d(d_in)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x):

        residual = x
        x = self.w_2(F.relu(self.w_1(x)))
        x = self.dropout(x)
        x += residual

        seq_len = x.size(1)
        x = self.batch_norm(x.view(-1, self.d_in)).view(-1, seq_len, self.d_in)
        return x

class EncoderLayer(nn.Module):
    ''' Compose with two layers '''

    def __init__(self, d_model, d_inner, n_head, d_k, d_v, dropout=0.1):
        super(EncoderLayer, self).__init__()
        self.slf_attn = MultiHeadAttention(n_head, d_model, d_k, d_v, dropout=dropout)
        self.pos_ffn = PositionwiseFeedForward(d_model, d_inner, dropout=dropout)

    def forward(self, enc_input, slf_attn_mask=None):
        enc_output, enc_slf_attn = self.slf_attn(
            enc_input, enc_input, enc_input, mask=slf_attn_mask)
        enc_output = self.pos_ffn(enc_output)
        return enc_output, enc_slf_attn

class Encoder(nn.Module):
    ''' A encoder model with self attention mechanism. '''

    def __init__(
            self, n_layers, n_head, d_k, d_v,
            d_model, d_inner, dropout=0., n_position=200):

        super().__init__()

        self.dropout = nn.Dropout(p=dropout)
        self.layer_stack = nn.ModuleList([
            EncoderLayer(d_model, d_inner, n_head, d_k, d_v, dropout=dropout)
            for _ in range(n_layers)])

    def forward(self, src_set, src_mask, return_attns=False):

        enc_slf_attn_list = []

        enc_output = src_set
        for enc_layer in self.layer_stack:
            enc_output, enc_slf_attn = enc_layer(enc_output, slf_attn_mask=src_mask)
            enc_slf_attn_list += [enc_slf_attn] if return_attns else []

        if return_attns:
            return enc_output, enc_slf_attn_list

        return enc_output

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

class Actor(pl.LightningModule):
    def __init__(self, degree, device=None):
        super().__init__()
        self.save_hyperparameters()
        self.degree = degree

        # embedder args
        self.d_input = 2
        self.d_model = 128
        self.embedder = Embedder(degree, self.d_input, self.d_model)

        # encoder args
        self.num_stacks = 3
        self.num_heads = 16
        self.d_k = 16
        self.d_v = 16
        # feedforward layer inner
        self.d_inner = 512

        self.encoder = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)

        # decoder args
        self.d_unit = 256
        self.d_query = 360
        self.conv1d_r = nn.Conv1d(self.d_model, self.d_unit, 1)
        self.conv1d_x = nn.Conv1d(self.d_model, self.d_unit, 1)
        self.conv1d_y = nn.Conv1d(self.d_model, self.d_unit, 1)

        self.start_ptr = Pointer(self.d_query, self.d_unit)
        self.q_l1 = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_l2 = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_l3 = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_lx = nn.Linear(self.d_model, self.d_query, bias=False)
        self.q_ly = nn.Linear(self.d_model, self.d_query, bias=False)
        self.relu = nn.ReLU()
        self.ctx_linear = nn.Linear(self.d_query, self.d_query, bias=False)

        self.ptr1 = Pointer(self.d_query, self.d_unit)
        self.ptr2 = Pointer(self.d_query, self.d_unit)

        self.train()

    def forward(self, inputs, deterministic=False):
        # encode encode encode
        inputs_tensor = inputs.detach().type(torch.float).to(self.device)
        embedings = self.embedder(inputs_tensor)
        encodings = self.encoder(embedings, None).permute(0, 2, 1)
        enc_r = self.conv1d_r(encodings).permute(0, 2, 1)
        enc_x = self.conv1d_x(encodings).permute(0, 2, 1)
        enc_y = self.conv1d_y(encodings).permute(0, 2, 1)
        encodings = encodings.permute(0, 2, 1)
        enc_xy = torch.cat([enc_x, enc_y], 1)

        batch_size = encodings.size()[0]

        visited = torch.zeros([batch_size, self.degree], dtype=torch.bool).to(self.device)

        indexes, log_probs = [], []

        # initial_idx = torch.zeros(batch_size, dtype=torch.long).to(self.device)
        start_logits = self.start_ptr(enc_r,
            torch.zeros([batch_size, self.d_query], dtype=torch.float).to(self.device), visited)
        distr = Categorical(logits=start_logits)
        if deterministic:
            _, start_idx = torch.max(start_logits, -1)
        else:
            start_idx = distr.sample()
        visited.scatter_(1, start_idx.unsqueeze(-1), True)
        log_probs.append(distr.log_prob(start_idx))

        q1 = encodings[torch.arange(batch_size), start_idx]
        q2 = q1
        qx = q1
        qy = q1

        context = torch.zeros([batch_size, self.d_query]).to(self.device)

        for step in range(self.degree - 1):
            residual = self.q_l1(q1) + self.q_l2(q2) + self.q_lx(qx) + self.q_ly(qy)
            context = torch.max(context, self.ctx_linear(self.relu(residual)))
            first_q = residual + context
            first_query = self.relu(first_q)
            # first_idx
            logits = self.ptr1(enc_r, first_query, visited)
            distr = Categorical(logits=logits)
            if deterministic:
                _, first_idx = torch.max(logits, -1)
            else:
                first_idx = distr.sample()
            log_probs.append(distr.log_prob(first_idx))

            # second_idx
            q3 = encodings[torch.arange(encodings.size(0)), first_idx]
            second_query = self.relu(first_q + self.q_l3(q3))

            unvisited = ~visited
            unvisited = torch.cat([unvisited, unvisited], -1)
            logits = self.ptr2(enc_xy, second_query, unvisited)
            distr = Categorical(logits=logits)
            if deterministic:
                _, idxs = torch.max(logits, -1)
            else:
                idxs = distr.sample()
            log_probs.append(distr.log_prob(idxs))

            second_idx = idxs % self.degree
            sec_dir = torch.div(idxs, self.degree, rounding_mode='floor')
            fir_dir = 1 - sec_dir

            with torch.no_grad():
                x_idx = first_idx * sec_dir + second_idx * fir_dir
                y_idx = first_idx * fir_dir + second_idx * sec_dir

            indexes.append(x_idx)
            indexes.append(y_idx)

            # update visited
            visited.scatter_(1, first_idx.unsqueeze(-1), True)

            # update query
            q1 = q3
            q2 = encodings[torch.arange(encodings.size(0)), second_idx]
            qx = encodings[torch.arange(encodings.size(0)), x_idx]
            qy = encodings[torch.arange(encodings.size(0)), y_idx]

        indexes = torch.stack(indexes, -1)
        log_probs = sum(log_probs)

        return indexes, log_probs

class Critic(pl.LightningModule):
    def __init__(self, degree, device=None):
        super().__init__()
        self.save_hyperparameters()
        self.degree = degree

        # embedder args
        self.d_input = 2
        self.d_model = 128

        # encoder args
        self.num_stacks = 3
        self.num_heads = 16
        self.d_k = 16
        self.d_v = 16
        self.d_inner = 512
        self.d_unit = 256

        self.crit_embedder = Embedder(degree, self.d_input, self.d_model)
        self.crit_encoder = Encoder(self.num_stacks, self.num_heads, self.d_k, self.d_v, self.d_model, self.d_inner)
        self.glimpse = Glimpse(self.d_model, self.d_unit)
        self.critic_l1 = nn.Linear(self.d_model, self.d_unit)
        self.critic_l2 = nn.Linear(self.d_unit, 1)
        self.relu = nn.ReLU()

        self.train()

    def forward(self, inputs, deterministic=False):
        inputs_tensor = torch.tensor(inputs, dtype=torch.float).to(self.device)
        critic_encode = self.crit_encoder(self.crit_embedder(inputs_tensor), None)
        glimpse = self.glimpse(critic_encode)
        critic_inner = self.relu(self.critic_l1(glimpse))
        predictions = self.relu(self.critic_l2(critic_inner)).squeeze(-1)

        return predictions