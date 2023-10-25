from inspect import isfunction
from functools import partial
import torch
import torch.nn as nn
import torch.nn.functional as F
import lightning.pytorch as pl
from einops import rearrange

def default(val, d):
    if val is not None:
        return val
    return d() if isfunction(d) else d
    
def zero_module(module):
    """
    Zero out the parameters of a module and return it.
    """
    for p in module.parameters():
        p.detach().zero_()
    return module

def Normalize(in_channels):
    return torch.nn.GroupNorm(num_groups=8, num_channels=in_channels, eps=1e-6, affine=True)
    
def split_at_index(dim, index, t):
    pre_slices = (slice(None),) * dim
    l = (*pre_slices, slice(None, index))
    r = (*pre_slices, slice(index, None))
    return t[l], t[r]
    
def causal_linear_attn(q, k, v, kv_mask = None, bucket_size = None, eps = 1e-3):
    b, h, n, e, dtype = *q.shape, q.dtype
    bucket_size = default(bucket_size, 64)
    bucket_size = max(bucket_size, 1)
    assert bucket_size == 0 or (n % bucket_size) == 0, f'sequence length {n} must be divisible by the bucket size {bucket_size} for causal linear attention'

    q = q.softmax(dim=-1)
    k = torch.exp(k).type(dtype).clone()

    q = q * e ** -0.5

    if kv_mask is not None:
        mask = kv_mask[:, None, :, None]
        k = k.masked_fill_(~mask, 0.)
        v = v.masked_fill_(~mask, 0.)
        del mask

    bucket_fn = lambda x: x.reshape(*x.shape[:-2], -1, bucket_size, e)
    b_q, b_k, b_v = map(bucket_fn, (q, k, v))

    b_k_sum = b_k.sum(dim=-2)
    b_k_cumsum = b_k_sum.cumsum(dim = -2).type(dtype)

    context = torch.einsum('bhund,bhune->bhude', b_k, b_v)
    context = context.cumsum(dim = -3).type(dtype)

    if bucket_size > 1:
        context = F.pad(context, (0, 0, 0, 0, 1, 0), value = 0.)
        context, _ = split_at_index(2, -1, context)

        b_k_cumsum = F.pad(b_k_cumsum, (0, 0, 1, 0), value = 0.)
        b_k_cumsum, _ = split_at_index(2, -1, b_k_cumsum)
    # print('b_q: ', b_q.shape)
    # print('b_k_cumsum: ', b_k_cumsum.shape)
    D_inv = 1. / torch.einsum('bhud,bhund->bhun', b_k_cumsum, b_q).clamp(min = eps)
    attn = torch.einsum('bhund,bhude,bhun->bhune', b_q, context, D_inv)
    return attn.reshape(*q.shape)
    
# feedforward
class GEGLU(nn.Module):
    def __init__(self, dim_in, dim_out):
        super().__init__()
        self.proj = nn.Linear(dim_in, dim_out * 2)

    def forward(self, x):
        x, gate = self.proj(x).chunk(2, dim=-1)
        return x * F.gelu(gate)

class FeedForward(nn.Module):
    def __init__(self, dim, dim_out=None, mult=4, glu=False, dropout=0.):
        super().__init__()
        inner_dim = int(dim * mult)
        dim_out = default(dim_out, dim)
        project_in = nn.Sequential(
            nn.Linear(dim, inner_dim),
            nn.GELU()
        ) if not glu else GEGLU(dim, inner_dim)

        self.net = nn.Sequential(
            project_in,
            nn.Dropout(dropout),
            nn.Linear(inner_dim, dim_out)
        )

    def forward(self, x):
        return self.net(x)

class SelfAttention(nn.Module):
    def __init__(self, query_dim, context_dim=None, heads=8, dim_head=64, dropout=0.):
        super().__init__()
        assert dim_head or (query_dim % heads) == 0, 'embedding dimension must be divisible by number of heads'
        dim_head = default(dim_head, query_dim // heads)
        context_dim = default(context_dim, query_dim)
        
        self.heads = heads
        self.dim_head = dim_head

        self.global_attn_heads = heads
        self.global_attn_fn = partial(causal_linear_attn, bucket_size = 1)

        self.to_q = nn.Linear(query_dim, dim_head * heads, bias = False)
        self.to_k = nn.Linear(context_dim, dim_head * heads, bias = False)
        self.to_v = nn.Linear(context_dim, dim_head * heads, bias = False)

        self.to_out = nn.Linear(dim_head * heads, query_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x, input_mask = None, context = None, context_mask = None):
        if context is None:
            q, k, v = (self.to_q(x), self.to_k(x), self.to_v(x))
        else:
            q, k, v = (self.to_q(x), self.to_k(context), self.to_v(context))
        b, t, e, h, dh = *q.shape, self.heads, self.dim_head

        merge_heads = lambda x: x.reshape(*x.shape[:2], -1, dh).transpose(1, 2)

        q, k, v = map(merge_heads, (q, k, v))

        out = []

        split_index_fn = partial(split_at_index, 1, 0)

        (lq, q), (lk, k), (lv, v) = map(split_index_fn, (q, k, v))

        has_local, has_global = map(lambda x: x.shape[1] > 0, (lq, q))

        if has_local:
            local_out = self.local_attn(lq, lk, lv, input_mask = input_mask)
            out.append(local_out)

        if has_global:
            kv_mask = context_mask
            global_out = self.global_attn_fn(q, k, v, kv_mask = kv_mask)
            out.append(global_out)

        attn = torch.cat(out, dim=1)
        attn = attn.transpose(1, 2).reshape(b, t, -1)
        return self.dropout(self.to_out(attn))

class BasicTransformerBlock(nn.Module):
    def __init__(self, dim, n_heads, d_head, dropout=0., context_dim=None, gated_ff=True):
        super().__init__()
        self.attn1 = SelfAttention(query_dim=dim, heads=n_heads, dim_head=d_head, dropout=dropout)  # is self-attention
        self.ff = FeedForward(dim, dropout=dropout, glu=gated_ff)
        self.attn2 = SelfAttention(query_dim=dim, context_dim=context_dim,
                                    heads=n_heads, dim_head=d_head, dropout=dropout)  # is cross-attn 
        self.norm1 = nn.LayerNorm(dim)
        self.norm2 = nn.LayerNorm(dim)
        self.norm3 = nn.LayerNorm(dim)

    def forward(self, x, context=None):
        x = self.attn1(self.norm1(x)) + x
        x = self.attn2(self.norm2(x), context=context) + x
        x = self.ff(self.norm3(x)) + x
        return x

class SpatialTransformer(nn.Module):
    """
    Transformer block for image-like data.
    First, project the input (aka embedding)
    and reshape to b, t, d.
    Then apply standard transformer action.
    Finally, reshape to image
    """
    def __init__(self, in_channels, n_heads=4, d_head=64,
                 depth=1, dropout=0., context_dim=None):
        super().__init__()
        inner_dim = n_heads * d_head
        self.norm = Normalize(in_channels)

        self.proj_in = nn.Conv2d(in_channels,
                                 inner_dim,
                                 kernel_size=1,
                                 stride=1,
                                 padding=0)

        self.transformer_blocks = nn.ModuleList(
            [BasicTransformerBlock(inner_dim, n_heads, d_head, dropout=dropout, context_dim=context_dim)
                for d in range(depth)]
        )

        self.proj_out = zero_module(nn.Conv2d(inner_dim,
                                              in_channels,
                                              kernel_size=1,
                                              stride=1,
                                              padding=0))

    def forward(self, x, context=None):
        # note: if no context is given, cross-attention defaults to self-attention
        b, c, h, w = x.shape
        x_in = x
        x = self.norm(x)
        x = self.proj_in(x)
        x = rearrange(x, 'b c h w -> b (h w) c')
        for block in self.transformer_blocks:
            x = block(x, context=context)
        x = rearrange(x, 'b (h w) c -> b c h w', h=h, w=w)
        x = self.proj_out(x)
        return x + x_in

class ResidualConvBlock(nn.Module):
    def __init__(
        self, in_channels: int, out_channels: int, z_channels: int, is_res: bool = False
    ) -> None:
        super().__init__()
        '''
        standard ResNet style convolutional block
        '''
        self.same_channels = in_channels==out_channels
        self.is_res = is_res
        self.conv1 = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, z_channels, 1, 1),
            nn.BatchNorm2d(out_channels),
            nn.GELU(),
        )
        self.conv2 = nn.Sequential(
            nn.Conv2d(out_channels, out_channels, z_channels, 1, 1),
            nn.BatchNorm2d(out_channels),
            nn.GELU(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.is_res:
            x1 = self.conv1(x)
            x2 = self.conv2(x1)
            # this adds on correct residual in case channels have increased
            if self.same_channels:
                out = x + x2
            else:
                out = x1 + x2 
            return out / 1.414
        else:
            x1 = self.conv1(x)
            x2 = self.conv2(x1)
            return x2

class UnetDown(nn.Module):
    def __init__(self, in_channels, out_channels, z_channels):
        super().__init__()
        '''
        process and downscale the image feature maps
        '''
        self.model = nn.Sequential(*[ResidualConvBlock(in_channels, out_channels, z_channels), nn.MaxPool2d(2)])

    def forward(self, x):
        return self.model(x)


class UnetUp(nn.Module):
    def __init__(self, in_channels, out_channels, z_channels):
        super().__init__()
        '''
        process and upscale the image feature maps
        '''
        layers = [
            nn.ConvTranspose2d(in_channels, out_channels, 2, 2),
            ResidualConvBlock(out_channels, out_channels, z_channels),
            ResidualConvBlock(out_channels, out_channels, z_channels),
        ]
        self.model = nn.Sequential(*layers)

    def forward(self, x, skip):
        x = torch.cat((x, skip), 1)
        x = self.model(x)
        return x


class EmbedFC(nn.Module):
    def __init__(self, input_dim, hidden_dim, emb_dim):
        super().__init__()
        '''
        generic one layer FC NN for embedding things  
        '''
        self.input_dim = input_dim
        layers = [
            nn.Linear(input_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, emb_dim),
        ]
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        x = x.view(-1, self.input_dim)
        return self.model(x)

class ResnetBlock(nn.Module):
    """Define a Resnet block"""

    def __init__(self, dim, padding_type, norm_layer, use_dropout, use_bias):
        """Initialize the Resnet block

        A resnet block is a conv block with skip connections
        We construct a conv block with build_conv_block function,
        and implement skip connections in <forward> function.
        Original Resnet paper: https://arxiv.org/pdf/1512.03385.pdf
        """
        super(ResnetBlock, self).__init__()
        self.conv_block = self.build_conv_block(dim, padding_type, norm_layer, use_dropout, use_bias)

    def build_conv_block(self, dim, padding_type, norm_layer, use_dropout, use_bias):
        """Construct a convolutional block.

        Parameters:
            dim (int)           -- the number of channels in the conv layer.
            padding_type (str)  -- the name of padding layer: reflect | replicate | zero
            norm_layer          -- normalization layer
            use_dropout (bool)  -- if use dropout layers.
            use_bias (bool)     -- if the conv layer uses bias or not

        Returns a conv block (with a conv layer, a normalization layer, and a non-linearity layer (ReLU))
        """
        conv_block = []
        p = 0
        if padding_type == 'reflect':
            conv_block += [nn.ReflectionPad2d(1)]
        elif padding_type == 'replicate':
            conv_block += [nn.ReplicationPad2d(1)]
        elif padding_type == 'zero':
            p = 1
        else:
            raise NotImplementedError('padding [%s] is not implemented' % padding_type)

        conv_block += [nn.Conv2d(dim, dim, kernel_size=3, padding=p, bias=use_bias), norm_layer(dim), nn.ReLU(True)]
        if use_dropout:
            conv_block += [nn.Dropout(0.25)]  # original dropout = 0.5

        p = 0
        if padding_type == 'reflect':
            conv_block += [nn.ReflectionPad2d(1)]
        elif padding_type == 'replicate':
            conv_block += [nn.ReplicationPad2d(1)]
        elif padding_type == 'zero':
            p = 1
        else:
            raise NotImplementedError('padding [%s] is not implemented' % padding_type)
        conv_block += [nn.Conv2d(dim, dim, kernel_size=3, padding=p, bias=use_bias), norm_layer(dim)]

        return nn.Sequential(*conv_block)

    def forward(self, x):
        """Forward function (with skip connections)"""
        out = x + self.conv_block(x)  # add skip connections
        return out

class ContextUnet(pl.LightningModule):
    def __init__(self, 
                 size, 
                 in_channels, 
                 c_in_channels, 
                 z_channels, 
                 n_feat=256, 
                 n_heads=4, 
                 dim_head=64, 
                 attention=True):
        super().__init__()
        self.save_hyperparameters()
        self.in_channels = in_channels
        self.n_feat = n_feat
        self.size = size
        self.attention = attention
        self.mask_value = torch.zeros([c_in_channels, size, size])

        self.init_conv = ResidualConvBlock(in_channels, n_feat, z_channels, is_res=True)

        self.down1 = UnetDown(n_feat, n_feat, z_channels)
        self.down2 = UnetDown(n_feat, 2*n_feat, z_channels)

        self.to_vec = nn.Sequential(nn.AvgPool2d(size//4), nn.GELU())

        self.timeembed1 = EmbedFC(1, 2*n_feat, 2*n_feat)
        self.timeembed2 = EmbedFC(1, 1*n_feat, 1*n_feat)
        if self.attention == True:
            self.conditionemb1 = EmbedFC(size*size*c_in_channels, 2*n_feat, (size//4)**2)
            self.conditionemb2 = EmbedFC(size*size*c_in_channels, 2*n_feat, (size//2)**2)
            self.st1_pc = SpatialTransformer(n_feat*2, n_heads=n_heads, d_head=dim_head, context_dim=1)
            self.st2_pc = SpatialTransformer(n_feat*1, n_heads=n_heads, d_head=dim_head, context_dim=1)
        else:
            self.conditionemb1 = EmbedFC(size*size*c_in_channels, 2*n_feat, 2*n_feat)
            self.conditionemb2 = EmbedFC(size*size*c_in_channels, 2*n_feat, 1*n_feat)

        self.up0 = nn.Sequential(
            nn.ConvTranspose2d(2*n_feat, 2*n_feat, size//4, size//4), # otherwise just have 2*n_feat
            nn.GroupNorm(8, 2*n_feat),
            nn.ReLU(),
        )
        
        self.up1 = UnetUp(2*2*n_feat, 1*n_feat, z_channels)
        self.up2 = UnetUp(2*1*n_feat, 1*n_feat, z_channels)
        self.out = nn.Sequential(
            nn.Conv2d(2*n_feat, n_feat, z_channels, 1, 1),
            nn.GroupNorm(8, n_feat),
            nn.ReLU(),
            nn.Conv2d(n_feat, self.in_channels, z_channels, 1, 1),
        )
        

    def forward(self, x, c, t, context_mask=None):
        # x is (noisy) connector image, c is capacity, p is pins, t is timestep 
        batch_size = x.shape[0]
        x = self.init_conv(x)
        down1 = self.down1(x)
        down2 = self.down2(down1)
        hiddenvec = self.to_vec(down2)

        # mask out context if context_mask == 1
        if context_mask is not None:
            for index in context_mask:
                c[index] *= self.mask_value.to(c.device)
                        
        temb1 = self.timeembed1(t).view(batch_size, 2*self.n_feat, 1, 1)
        temb2 = self.timeembed2(t).view(batch_size, 1*self.n_feat, 1, 1)
        
        if self.attention == True:
            cemb1 = self.conditionemb1(c).view(batch_size, (self.size//4)**2, 1)
            cemb2 = self.conditionemb2(c).view(batch_size, (self.size//2)**2, 1)
        else:
            cemb1 = self.conditionemb1(c).view(-1, self.n_feat * 2, 1, 1)
            cemb2 = self.conditionemb2(c).view(-1, self.n_feat, 1, 1)
                
        up1 = self.up0(hiddenvec)
        
        if self.attention == True:
            up1_ca_pc = self.st1_pc(up1, cemb1)
            up2 = self.up1(up1_ca_pc + temb1, down2)
            up2_ca_pc = self.st2_pc(up2, cemb2)
            up3 = self.up2(up2_ca_pc + temb2, down1)
        else:
            up2 = self.up1(up1*cemb1 + temb1, down2)
            up3 = self.up2(up2*cemb2 + temb2, down1)
        
        out = self.out(torch.cat((up3, x), 1))
        
        return out
