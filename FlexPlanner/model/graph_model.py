import torch
from torch import nn
from einops import rearrange, repeat
from tianshou.data import Batch


class PositionalEncoding(nn.Module):
    def __init__(self, max_len:int, dim:int, trainable:bool=True):
        super().__init__()
        if trainable:
            self.pos_enc = nn.Parameter(torch.randn(max_len, dim))
        else:
            pos_enc = torch.zeros(max_len, dim)
            pos = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
            div_term = torch.exp(torch.arange(0, dim, 2).float() * (-torch.log(torch.tensor(10000.0)) / dim))
            pos_enc[:, 0::2] = torch.sin(pos * div_term)
            pos_enc[:, 1::2] = torch.cos(pos * div_term)
            self.register_buffer('pos_enc', pos_enc)

    def forward(self, x:torch.Tensor, order:torch.LongTensor) -> torch.Tensor:
        """order: [B,L]"""
        if order is None:
            return x
        
        batch_size, seq_len = order.shape
        pe = self.pos_enc[order.flatten()]
        pe = pe.view(batch_size, seq_len, -1)
        return pe + x


class Transformer(nn.Module):
    """
    This is the backbone to process graph feature in the shared encoder.
    """
    def __init__(self, in_dim:int, hid_dim:int, out_dim:int, n_layers:int, episode_len:int):
        super().__init__()
        self.fc0 = nn.Linear(in_dim, hid_dim)
        self.fc1 = nn.Linear(hid_dim, out_dim)
        self.fc2 = nn.Linear(hid_dim, out_dim)
        self.nhead = 4

        enc_layer = nn.TransformerEncoderLayer(d_model=hid_dim, nhead=self.nhead, dim_feedforward=hid_dim*2, dropout=0.0, batch_first=True)
        self.transformer = nn.TransformerEncoder(enc_layer, num_layers=n_layers)
        self.invalid_node = nn.Parameter(torch.randn(1, hid_dim))

        self.pe = PositionalEncoding(episode_len+1, hid_dim, trainable=False)

    def forward(self, x:torch.Tensor, adj_mat:torch.IntTensor, idx:torch.LongTensor, order:torch.LongTensor) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Input:
            x: the input tensor, is a torch.Tensor, shape = [B, L, D].
        Output:
            output: global and local embeddings, both are torch.Tensor, shape = [B, D].
        """
        x = self.fc0(x)
        x = self.pe(x, order)
        mask = (1-adj_mat).bool()
        mask = repeat(mask, 'b i j -> b n i j', n=self.nhead)
        mask = rearrange(mask, 'b n i j -> (b n) i j')
        x = self.transformer(x, mask)

        if (idx == -1).sum() > 0:
            pass
        
        emb_local = x[torch.arange(x.size(0)), idx]
        emb_local[idx == -1] = self.invalid_node
        emb_global = x.mean(dim=1)

        emb_local = self.fc1(emb_local)
        emb_global = self.fc2(emb_global)
        return emb_global, emb_local


def create_graph_model(graph:int, shared_encoder_output_dim:int, episode_len:int, hid_dim:int=32, transformer_out_dim:int=32) -> tuple[Transformer, nn.Linear]:
    """
    Arguments:
        @graph:int.
            0: no graph.
            1: use transformer (GNN) to extract global and local graph features, then add them to the output.
            2: use transformer (GNN) to extract global and local graph features, then concatenate them to the output.
            3: use transformer (GNN) to extract global and local graph features, then concatenate (global+local) to the output.
    """
    assert graph > 0

    in_dim = 7

    if graph == 1:
        print("[WARNING] graph = {} is not recommended.".format(graph))
        transformer_out_dim = shared_encoder_output_dim
        fc_graph = nn.Identity()

    elif graph == 2:
        fc_graph = nn.Sequential(
            nn.Linear(transformer_out_dim*2+shared_encoder_output_dim, shared_encoder_output_dim),
            nn.ReLU(),
            nn.Linear(shared_encoder_output_dim, shared_encoder_output_dim),
        )

    elif graph == 3:
        print("[WARNING] graph = {} is not recommended.".format(graph))
        fc_graph = nn.Sequential(
            nn.Linear(transformer_out_dim+shared_encoder_output_dim, shared_encoder_output_dim),
            nn.ReLU(),
            nn.Linear(shared_encoder_output_dim, shared_encoder_output_dim),
        )

    else:
        raise ValueError(f"graph should be 0, 1, 2, 3, got {graph}")
    
    transformer = Transformer(in_dim, hid_dim, transformer_out_dim, n_layers=2, episode_len=episode_len)
    return transformer, fc_graph


def forward_graph_model(graph:int, graph_model:Transformer, fc_graph:nn.Module, x:torch.Tensor, graph_data_batch:Batch) -> torch.Tensor:
    """
    Arguments:
        @graph:int.
            0: no graph.
            1: use transformer (GNN) to extract global and local graph features, then add them to the output.
            2: use transformer (GNN) to extract global and local graph features, then concatenate them to the output.
            3: use transformer (GNN) to extract global and local graph features, then concatenate (global+local) to the output.
    """
    assert graph > 0
    graph_data = [graph_data_batch[k] for k in ['x', 'y', 'z', 'w', 'h', 'area', 'placed']]
    graph_data = torch.stack(graph_data, dim=-1).to(x.device)
    adj_mat = graph_data_batch['adj_mat_mov'].to(x.device)
    curr_node = graph_data_batch['idx'].long()
    order = graph_data_batch['order'].long() if 'order' in graph_data_batch.keys() else None

    emb_global, emb_local = graph_model(graph_data, adj_mat, curr_node, order)

    if graph == 1:
        y = x + emb_global + emb_local

    elif graph == 2:
        y = fc_graph(torch.cat([x, emb_global, emb_local], dim=-1))

    elif graph == 3:
        y = fc_graph(torch.cat([x, emb_local+emb_global], dim=-1))

    return y