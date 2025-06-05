# --------------------------------------------------------
# ref:
# https://github.com/microsoft/unilm/tree/master/beit
# https://github.com/rwightman/pytorch-image-models/tree/master/timm
# https://github.com/facebookresearch/deit/
# https://github.com/facebookresearch/dino
# --------------------------------------------------------'
import math
from functools import partial

import torch
import torch.nn as nn
import torch.nn.functional as F
from timm.models.layers import drop_path, to_2tuple, trunc_normal_
from timm.models.registry import register_model
from copy import deepcopy


def _cfg(url='', **kwargs):
    return {
        'url': url,
        'num_classes': 1000, 'input_size': (3, 224, 224), 'pool_size': None,
        'crop_pct': .9, 'interpolation': 'bicubic',
        'mean': (0.5, 0.5, 0.5), 'std': (0.5, 0.5, 0.5),
        **kwargs
    }


def convert_module_to_f32(l):
    """
    Convert primitive modules to float32, undoing convert_module_to_f16().
    """
    if isinstance(l, (nn.Conv1d, nn.Conv2d, nn.Conv3d, nn.Linear)):
        l.weight.data = l.weight.data.float()
        if l.bias is not None:
            l.bias.data = l.bias.data.float()


def convert_module_to_f16(l):
    """
    Convert primitive modules to float16.
    """
    if isinstance(l, (nn.Conv1d, nn.Conv2d, nn.Conv3d, nn.Linear)):
        l.weight.data = l.weight.data.half()
        if l.bias is not None:
            l.bias.data = l.bias.data.half()


class DropPath(nn.Module):
    """Drop paths (Stochastic Depth) per sample  (when applied in main path of residual blocks).
    """
    def __init__(self, drop_prob=None):
        super(DropPath, self).__init__()
        self.drop_prob = drop_prob

    def forward(self, x):
        return drop_path(x, self.drop_prob, self.training)
    
    def extra_repr(self) -> str:
        return 'p={}'.format(self.drop_prob)


class Mlp(nn.Module):
    def __init__(self, in_features, hidden_features=None, out_features=None, act_layer=nn.GELU, drop=0.):
        super().__init__()
        out_features = out_features or in_features
        hidden_features = hidden_features or in_features
        self.fc1 = nn.Linear(in_features, hidden_features)
        self.act = act_layer()
        self.fc2 = nn.Linear(hidden_features, out_features)
        self.drop = nn.Dropout(drop)

    def forward(self, x):
        x = self.fc1(x)
        x = self.act(x)
        # x = self.drop(x)
        # commit this for the orignal BERT implement 
        x = self.fc2(x)
        x = self.drop(x)
        return x


class Attention(nn.Module):
    def __init__(
            self, dim, num_heads=8, qkv_bias=False, qk_scale=None, attn_drop=0.,
            proj_drop=0., window_size=None, attn_head_dim=None, num_extra_tokens=1):
        super().__init__()
        self.num_heads = num_heads
        head_dim = dim // num_heads
        if attn_head_dim is not None:
            head_dim = attn_head_dim
        all_head_dim = head_dim * self.num_heads
        self.scale = qk_scale or head_dim ** -0.5

        self.qkv = nn.Linear(dim, all_head_dim * 3, bias=False)
        if qkv_bias:
            self.q_bias = nn.Parameter(torch.zeros(all_head_dim))
            self.v_bias = nn.Parameter(torch.zeros(all_head_dim))
        else:
            self.q_bias = None
            self.v_bias = None

        if window_size:
            self.window_size = window_size

            # self.register_buffer("relative_position_index", relative_position_index)
            assert num_extra_tokens in [0, 1, 2]
            extra_relative_distance = num_extra_tokens * (num_extra_tokens + 2)
            self.num_extra_tokens = num_extra_tokens

            self.num_relative_distance = (2 * window_size[0] - 1) * (2 * window_size[1] - 1) \
                + extra_relative_distance
            self.relative_position_bias_table = nn.Parameter(
                torch.zeros(self.num_relative_distance, num_heads))  # 2*Wh-1 * 2*Ww-1, nH

            # get pair-wise relative position index for each token inside the window
            coords_h = torch.arange(window_size[0])
            coords_w = torch.arange(window_size[1])
            coords = torch.stack(torch.meshgrid([coords_h, coords_w]))  # 2, Wh, Ww
            coords_flatten = torch.flatten(coords, 1)  # 2, Wh*Ww
            relative_coords = coords_flatten[:, :, None] - coords_flatten[:, None, :]  # 2, Wh*Ww, Wh*Ww
            relative_coords = relative_coords.permute(1, 2, 0).contiguous()  # Wh*Ww, Wh*Ww, 2
            relative_coords[:, :, 0] += window_size[0] - 1  # shift to start from 0
            relative_coords[:, :, 1] += window_size[1] - 1
            relative_coords[:, :, 0] *= 2 * window_size[1] - 1
            relative_position_index = \
                torch.zeros(size=(window_size[0] * window_size[1] + num_extra_tokens,) * 2, dtype=relative_coords.dtype)
            relative_position_index[num_extra_tokens:, num_extra_tokens:] = relative_coords.sum(-1)  # Wh*Ww, Wh*Ww
            if num_extra_tokens == 1:
                relative_position_index[0, 0:] = self.num_relative_distance - 3
                relative_position_index[0:, 0] = self.num_relative_distance - 2
                relative_position_index[0, 0]  = self.num_relative_distance - 1
            elif num_extra_tokens == 2:
                relative_position_index[1, 1]  = self.num_relative_distance - 8
                relative_position_index[1, 0]  = self.num_relative_distance - 7
                relative_position_index[0, 1]  = self.num_relative_distance - 6
                relative_position_index[0, 2:] = self.num_relative_distance - 5
                relative_position_index[2:, 0] = self.num_relative_distance - 4
                relative_position_index[1, 2:] = self.num_relative_distance - 3
                relative_position_index[2:, 1] = self.num_relative_distance - 2
                relative_position_index[0, 0]  = self.num_relative_distance - 1

            self.register_buffer("relative_position_index", relative_position_index)
        else:
            self.window_size = None
            self.relative_position_bias_table = None
            self.relative_position_index = None

        self.attn_drop = nn.Dropout(attn_drop)
        self.proj = nn.Linear(all_head_dim, dim)
        self.proj_drop = nn.Dropout(proj_drop)

    def forward(self, x, rel_pos_bias=None):
        B, N, C = x.shape
        qkv_bias = None
        if self.q_bias is not None:
            qkv_bias = torch.cat((self.q_bias, torch.zeros_like(self.v_bias, requires_grad=False), self.v_bias))
        # qkv = self.qkv(x).reshape(B, N, 3, self.num_heads, C // self.num_heads).permute(2, 0, 3, 1, 4)
        qkv = F.linear(input=x, weight=self.qkv.weight, bias=qkv_bias.type(x.dtype))
        qkv = qkv.reshape(B, N, 3, self.num_heads, -1).permute(2, 0, 3, 1, 4)
        q, k, v = qkv[0], qkv[1], qkv[2]   # make torchscript happy (cannot use tensor as tuple)

        q = q * self.scale
        # attn = (q @ k.transpose(-2, -1))
        attn = (q.float() @ k.float().transpose(-2, -1))

        if self.relative_position_bias_table is not None:
            relative_position_bias = \
                self.relative_position_bias_table[self.relative_position_index.view(-1)].view(
                    self.window_size[0] * self.window_size[1] + self.num_extra_tokens,
                    self.window_size[0] * self.window_size[1] + self.num_extra_tokens, -1)  # Wh*Ww,Wh*Ww,nH
            relative_position_bias = relative_position_bias.permute(2, 0, 1).contiguous()  # nH, Wh*Ww, Wh*Ww
            attn = attn + relative_position_bias.unsqueeze(0)

        if rel_pos_bias is not None:
            attn = attn + rel_pos_bias
        
        attn = attn.softmax(dim=-1)
        attn = self.attn_drop(attn)

        x = (attn @ v.type(attn.dtype)).transpose(1, 2).reshape(B, N, -1).type(x.dtype)
        x = self.proj(x)
        x = self.proj_drop(x)
        return x


class Block(nn.Module):

    def __init__(self, dim, num_heads, mlp_ratio=4., qkv_bias=False, qk_scale=None, drop=0., attn_drop=0.,
                 drop_path=0., init_values=None, act_layer=nn.GELU, norm_layer=nn.LayerNorm,
                 window_size=None, attn_head_dim=None, num_extra_tokens=1):
        super().__init__()
        self.norm1 = norm_layer(dim)
        self.attn = Attention(
            dim, num_heads=num_heads, qkv_bias=qkv_bias, qk_scale=qk_scale, num_extra_tokens=num_extra_tokens,
            attn_drop=attn_drop, proj_drop=drop, window_size=window_size, attn_head_dim=attn_head_dim)
        # NOTE: drop path for stochastic depth, we shall see if this is better than dropout here
        self.drop_path = DropPath(drop_path) if drop_path > 0. else nn.Identity()
        self.norm2 = norm_layer(dim)
        mlp_hidden_dim = int(dim * mlp_ratio)
        self.mlp = Mlp(in_features=dim, hidden_features=mlp_hidden_dim, act_layer=act_layer, drop=drop)

        if init_values is not None and init_values > 0:
            self.gamma_1 = nn.Parameter(init_values * torch.ones((dim)),requires_grad=True)
            self.gamma_2 = nn.Parameter(init_values * torch.ones((dim)),requires_grad=True)
        else:
            self.gamma_1, self.gamma_2 = None, None

    def forward(self, x, rel_pos_bias=None):
        if self.gamma_1 is None:
            x = x + self.drop_path(self.attn(self.norm1(x.float()).type(x.dtype), rel_pos_bias=rel_pos_bias))
            x = x + self.drop_path(self.mlp(self.norm2(x.float()).type(x.dtype)))
        else:
            x = x + self.drop_path(self.gamma_1 * self.attn(self.norm1(x.float()).type(x.dtype), rel_pos_bias=rel_pos_bias))
            x = x + self.drop_path(self.gamma_2 * self.mlp(self.norm2(x.float()).type(x.dtype)))
        return x


class PatchEmbed(nn.Module):
    """ Image to Patch Embedding
    """
    def __init__(self, img_size=224, patch_size=16, in_chans=3, embed_dim=768):
        super().__init__()
        img_size = to_2tuple(img_size)
        patch_size = to_2tuple(patch_size)
        num_patches = (img_size[1] // patch_size[1]) * (img_size[0] // patch_size[0])
        self.patch_shape = (img_size[0] // patch_size[0], img_size[1] // patch_size[1])
        self.img_size = img_size
        self.patch_size = patch_size
        self.num_patches = num_patches

        self.proj = nn.Conv2d(in_chans, embed_dim, kernel_size=patch_size, stride=patch_size)

    def forward(self, x, **kwargs):
        B, C, H, W = x.shape
        # FIXME look at relaxing size constraints
        assert H == self.img_size[0] and W == self.img_size[1], \
            f"Input image size ({H}*{W}) doesn't match model ({self.img_size[0]}*{self.img_size[1]})."
        x = self.proj(x).flatten(2).transpose(1, 2)
        return x


class RelativePositionBias(nn.Module):

    def __init__(self, window_size, num_heads, num_extra_tokens=1):
        super().__init__()
        self.window_size = window_size

        # cls to token & token 2 cls & cls to cls
        # if having an extra time token
        # time to token & token 2 time & time to time & cls to time & time to cls
        assert num_extra_tokens in [0, 1, 2]
        extra_relative_distance = num_extra_tokens * (num_extra_tokens + 2)
        self.num_extra_tokens = num_extra_tokens

        self.num_relative_distance = (2 * window_size[0] - 1) * (2 * window_size[1] - 1) \
            + extra_relative_distance
        self.relative_position_bias_table = nn.Parameter(
            torch.zeros(self.num_relative_distance, num_heads))  # 2*Wh-1 * 2*Ww-1, nH

        # get pair-wise relative position index for each token inside the window
        coords_h = torch.arange(window_size[0])
        coords_w = torch.arange(window_size[1])
        coords = torch.stack(torch.meshgrid([coords_h, coords_w]))  # 2, Wh, Ww
        coords_flatten = torch.flatten(coords, 1)  # 2, Wh*Ww
        relative_coords = coords_flatten[:, :, None] - coords_flatten[:, None, :]  # 2, Wh*Ww, Wh*Ww
        relative_coords = relative_coords.permute(1, 2, 0).contiguous()  # Wh*Ww, Wh*Ww, 2
        relative_coords[:, :, 0] += window_size[0] - 1  # shift to start from 0
        relative_coords[:, :, 1] += window_size[1] - 1
        relative_coords[:, :, 0] *= 2 * window_size[1] - 1
        relative_position_index = \
            torch.zeros(size=(window_size[0] * window_size[1] + num_extra_tokens,) * 2, dtype=relative_coords.dtype)
        relative_position_index[num_extra_tokens:, num_extra_tokens:] = relative_coords.sum(-1)  # Wh*Ww, Wh*Ww
        if num_extra_tokens == 1:
            relative_position_index[0, 0:] = self.num_relative_distance - 3
            relative_position_index[0:, 0] = self.num_relative_distance - 2
            relative_position_index[0, 0]  = self.num_relative_distance - 1
        elif num_extra_tokens == 2:
            relative_position_index[1, 1]  = self.num_relative_distance - 8
            relative_position_index[1, 0]  = self.num_relative_distance - 7
            relative_position_index[0, 1]  = self.num_relative_distance - 6
            relative_position_index[0, 2:] = self.num_relative_distance - 5
            relative_position_index[2:, 0] = self.num_relative_distance - 4
            relative_position_index[1, 2:] = self.num_relative_distance - 3
            relative_position_index[2:, 1] = self.num_relative_distance - 2
            relative_position_index[0, 0]  = self.num_relative_distance - 1

        self.register_buffer("relative_position_index", relative_position_index)

        # trunc_normal_(self.relative_position_bias_table, std=.02)

    def forward(self):
        relative_position_bias = \
            self.relative_position_bias_table[self.relative_position_index.view(-1)].view(
                self.window_size[0] * self.window_size[1] + self.num_extra_tokens,
                self.window_size[0] * self.window_size[1] + self.num_extra_tokens, -1)  # Wh*Ww,Wh*Ww,nH
        return relative_position_bias.permute(2, 0, 1).contiguous()  # nH, Wh*Ww, Wh*Ww


class VisionTransformer(nn.Module):
    """ Vision Transformer with support for patch or hybrid CNN input stage
    """
    def __init__(self, img_size=224, patch_size=16, in_chans=3, num_classes=1000, embed_dim=768, depth=12,
                 num_heads=12, mlp_ratio=4., qkv_bias=False, qk_scale=None, drop_rate=0., attn_drop_rate=0.,
                 drop_path_rate=0., norm_layer=nn.LayerNorm, init_values=None,
                 use_abs_pos_emb=True, use_rel_pos_bias=False, use_shared_rel_pos_bias=False,
                 use_mean_pooling=True, init_scale=0.001, use_conv_last=False, num_steps=4000,
                 learn_sigma=False, use_fp16=True, drop_label_prob=0.0, **kwargs):
        super().__init__()
        self.num_classes = num_classes
        self.num_features = self.embed_dim = embed_dim  # num_features for consistency with other models
        self.patch_size = patch_size
        self.img_size = img_size
        self.use_conv_last = use_conv_last
        self.drop_label_prob = drop_label_prob

        self.patch_embed = PatchEmbed(
            img_size=img_size, patch_size=patch_size, in_chans=in_chans, embed_dim=embed_dim)
        num_patches = self.patch_embed.num_patches

        self.num_extra_tokens = 1
        if num_classes > 0:
            self.num_extra_tokens += 1
        if use_abs_pos_emb:
            self.pos_embed = nn.Parameter(torch.zeros(1, num_patches + self.num_extra_tokens, embed_dim))
        else:
            self.pos_embed = None
        self.pos_drop = nn.Dropout(p=drop_rate)

        if use_shared_rel_pos_bias:
            self.rel_pos_bias = RelativePositionBias(
                window_size=self.patch_embed.patch_shape, 
                num_heads=num_heads, num_extra_tokens=self.num_extra_tokens)
        else:
            self.rel_pos_bias = None

        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, depth)]  # stochastic depth decay rule
        self.use_rel_pos_bias = use_rel_pos_bias
        self.blocks = nn.ModuleList([
            Block(
                dim=embed_dim, num_heads=num_heads, mlp_ratio=mlp_ratio, qkv_bias=qkv_bias, qk_scale=qk_scale,
                drop=drop_rate, attn_drop=attn_drop_rate, drop_path=dpr[i], norm_layer=norm_layer,
                init_values=init_values, num_extra_tokens=self.num_extra_tokens,
                window_size=self.patch_embed.patch_shape if use_rel_pos_bias else None)
            for i in range(depth)])
        self.norm = nn.Identity() if use_mean_pooling else norm_layer(embed_dim)

        if learn_sigma:
            self.out_dim = in_chans * 2
        else:
            self.out_dim = in_chans
        
        output_channels = self.out_dim * patch_size ** 2
        self.linear_projection = nn.Linear(embed_dim, output_channels)
        self.to_pixel = nn.Conv2d(self.out_dim, self.out_dim, kernel_size=3, stride=1, padding=1) if self.use_conv_last else None

        if self.pos_embed is not None:
            trunc_normal_(self.pos_embed, std=.02)

        if isinstance(self.linear_projection, nn.Linear):
            trunc_normal_(self.linear_projection.weight, std=.02)
        self.apply(self._init_weights)
        self.fix_init_weight()

        if isinstance(self.linear_projection, nn.Linear):
            self.linear_projection.weight.data.mul_(init_scale)
            self.linear_projection.bias.data.mul_(init_scale)

        # class embedding
        if num_classes > 0:
            self.class_embedding = nn.Embedding(
                num_classes + int(self.drop_label_prob > 0), embed_dim)
        else:
            self.class_embedding = None
        
        # timestep embedding
        self.time_embedding = nn.Embedding(num_steps, embed_dim)

        self.dtype = torch.float16 if use_fp16 else torch.float32

    def fix_init_weight(self):
        def rescale(param, layer_id):
            param.div_(math.sqrt(2.0 * layer_id))

        for layer_id, layer in enumerate(self.blocks):
            rescale(layer.attn.proj.weight.data, layer_id + 1)
            rescale(layer.mlp.fc2.weight.data, layer_id + 1)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            trunc_normal_(m.weight, std=.02)
            if isinstance(m, nn.Linear) and m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)

    def get_num_layers(self):
        return len(self.blocks)

    @torch.jit.ignore
    def no_weight_decay(self):
        return {'pos_embed', 'cls_token'}

    def get_classifier(self):
        return self.head

    def reset_classifier(self, num_classes, global_pool=''):
        self.num_classes = num_classes
        self.head = nn.Linear(self.embed_dim, num_classes) if num_classes > 0 else nn.Identity()

    def convert_to_fp16(self):
        """
        Convert the torso of the model to float16.
        """
        self.apply(convert_module_to_f16)
        self.patch_embed.apply(convert_module_to_f32)
        self.linear_projection.apply(convert_module_to_f32)
        if self.to_pixel is not None:
            self.to_pixel.apply(convert_module_to_f32)

        for blk in self.blocks:
            if blk.attn.relative_position_bias_table is not None:
                blk.attn.relative_position_bias_table = blk.attn.relative_position_bias_table.float()


    def convert_to_fp32(self):
        """
        Convert the torso of the model to float32.
        """
        self.apply(convert_module_to_f32)
    
    def token_drop(self, labels, force_drop_ids=None):
        """
        Drops labels to enable classifier-free guidance.
        """
        if force_drop_ids is None:
            drop_ids = torch.rand(labels.shape[0]).to(labels.device) < self.drop_label_prob
        else:
            drop_ids = force_drop_ids == 1
        labels = torch.where(drop_ids, self.num_classes, labels)
        return labels
    
    def forward_features(self, x, timesteps, y=None, force_drop_ids=None):
        
        x = self.patch_embed(x)
        B, L, C = x.size()

        time_tokens = self.time_embedding(timesteps).unsqueeze(1)

        if y is not None:
            use_dropout = self.drop_label_prob > 0 and self.training
            if use_dropout or (force_drop_ids is not None):
                y = self.token_drop(y, force_drop_ids)
            cls_tokens = self.class_embedding(y).unsqueeze(1)
            x = torch.cat((time_tokens, cls_tokens, x), dim=1)
        else:
            x = torch.cat((time_tokens, x), dim=1)

        if self.pos_embed is not None:
            x = x + self.pos_embed
        x = self.pos_drop(x)

        rel_pos_bias = self.rel_pos_bias() if self.rel_pos_bias is not None else None
        x = x.type(self.dtype)
        for blk in self.blocks:
            x = blk(x, rel_pos_bias=rel_pos_bias)

        x = self.norm(x)
        return x

    def forward(self, x, timesteps, y=None, force_drop_ids=None, **kwargs):
        input_dtype = x.dtype
        x = self.forward_features(x, timesteps, y=y, force_drop_ids=force_drop_ids)
        B, L, C = x.shape

        x = x.type(input_dtype)
        x = self.linear_projection(x[:, self.num_extra_tokens:, :])

        p = self.patch_size
        h = w = int(x.shape[1]**.5)
        assert h * w == x.shape[1]
        
        x = x.reshape(shape=(x.shape[0], h, w, p, p, self.out_dim))
        x = torch.einsum('nhwpqc->nchpwq', x)
        imgs = x.reshape(shape=(x.shape[0], self.out_dim, h * p, h * p))

        if self.to_pixel is not None:
            imgs = self.to_pixel(imgs)

        return imgs

    def forward_with_cfg(self, x, timesteps, y=None, classifier_free_scale=1.0):
        input_dtype = x.dtype
        # import pdb; pdb.set_trace()
        half = x[: len(x) // 2]
        combined = torch.cat((half, half), dim=0)
        x = self.forward_features(combined, timesteps, y=y, force_drop_ids=None)
        B, L, C = x.shape
        # x = self.head(x)
        x = x.type(input_dtype)
        x = self.linear_projection(x[:, self.num_extra_tokens:, :])

        p = self.patch_size
        h = w = int(x.shape[1]**.5)
        assert h * w == x.shape[1]
        
        x = x.reshape(shape=(x.shape[0], h, w, p, p, self.out_dim))
        x = torch.einsum('nhwpqc->nchpwq', x)
        imgs = x.reshape(shape=(x.shape[0], self.out_dim, h * p, h * p))

        if self.to_pixel is not None:
            imgs = self.to_pixel(imgs)

        _cond, _uncond = torch.split(imgs, len(imgs) // 2, dim=0)
        _half = _uncond + classifier_free_scale * (_cond - _uncond)

        imgs = torch.cat((_half, _half), dim=0)

        return imgs

    def get_intermediate_layers(self, x, timesteps, y=None):
        x = self.patch_embed(x)
        batch_size, seq_len, _ = x.size()

        # cls_tokens = self.cls_token.expand(batch_size, -1, -1)  # stole cls_tokens impl from Phil Wang, thanks
        # x = torch.cat((cls_tokens, x), dim=1)
        cls_tokens = self.class_embedding(y).unsqueeze(1)
        time_tokens = self.time_embedding(timesteps).unsqueeze(1)

        x = torch.cat((time_tokens, cls_tokens, x), dim=1)        
        
        if self.pos_embed is not None:
            x = x + self.pos_embed
        x = self.pos_drop(x)

        features = []
        rel_pos_bias = self.rel_pos_bias() if self.rel_pos_bias is not None else None
        for blk in self.blocks:
            x = blk(x, rel_pos_bias)
            features.append(x)

        return features


@register_model
def vit_base_patch2_32(pretrained=False, **kwargs):
    new_kwargs = kwargs
    for key in ['img_size', 'patch_size', 'embed_dim', 'num_heads', 'mlp_ratio', 'qkv_bias', 'norm_layer']:
        if key in new_kwargs:
            new_kwargs.pop(key)
    if not kwargs['class_cond']:
        kwargs['num_classes'] = -1
    model = VisionTransformer(
        img_size=32,
        patch_size=2, embed_dim=768, num_heads=12, mlp_ratio=4, qkv_bias=True,
        norm_layer=partial(nn.LayerNorm, eps=1e-6), **new_kwargs)
    model.default_cfg = _cfg()
    return model


@register_model
def vit_large_patch2_32(pretrained=False, **kwargs):
    new_kwargs = deepcopy(kwargs)
    for key in ['img_size', 'patch_size', 'embed_dim', 'num_heads', 'mlp_ratio', 'qkv_bias', 'norm_layer']:
        if key in new_kwargs:
            new_kwargs.pop(key)
    if not kwargs['class_cond']:
        new_kwargs['num_classes'] = -1
    model = VisionTransformer(
        img_size=32,
        patch_size=2, embed_dim=1024, num_heads=16, mlp_ratio=4, qkv_bias=True,
        norm_layer=partial(nn.LayerNorm, eps=1e-6), **new_kwargs)
    model.default_cfg = _cfg()
    return model


@register_model
def vit_large_patch4_64(pretrained=False, **kwargs):
    new_kwargs = deepcopy(kwargs)
    for key in ['img_size', 'patch_size', 'embed_dim', 'num_heads', 'mlp_ratio', 'qkv_bias', 'norm_layer']:
        if key in new_kwargs:
            new_kwargs.pop(key)
    if not kwargs['class_cond']:
        new_kwargs['num_classes'] = -1
    model = VisionTransformer(
        img_size=64,
        patch_size=4, embed_dim=1024, num_heads=16, mlp_ratio=4, qkv_bias=True,
        norm_layer=partial(nn.LayerNorm, eps=1e-6), **new_kwargs)
    model.default_cfg = _cfg()
    return model


@register_model
def vit_xl_patch2_32(pretrained=False, **kwargs):
    new_kwargs = deepcopy(kwargs)
    for key in ['img_size', 'patch_size', 'embed_dim', 'num_heads', 'mlp_ratio', 'qkv_bias', 'norm_layer']:
        if key in new_kwargs:
            new_kwargs.pop(key)
    if not kwargs['class_cond']:
        new_kwargs['num_classes'] = -1
    model = VisionTransformer(
        img_size=32,
        patch_size=2, embed_dim=1152, num_heads=16, mlp_ratio=4, qkv_bias=True,
        norm_layer=partial(nn.LayerNorm, eps=1e-6), **new_kwargs)
    model.default_cfg = _cfg()
    return model


@register_model
def vit_base_arbitrary(pretrained=False, **kwargs):
    new_kwargs = kwargs
    assert 'img_size' in new_kwargs and 'patch_size' in new_kwargs
    assert new_kwargs['img_size'] % new_kwargs['patch_size'] == 0

    for key in ['embed_dim', 'num_heads', 'mlp_ratio', 'qkv_bias', 'norm_layer']:
        if key in new_kwargs:
            new_kwargs.pop(key)
    if not new_kwargs['class_cond']:
        new_kwargs['num_classes'] = -1
    model = VisionTransformer(
        img_size=new_kwargs['img_size'],
        patch_size=new_kwargs['patch_size'], 
        embed_dim=768, num_heads=12, mlp_ratio=4, qkv_bias=True,
        norm_layer=partial(nn.LayerNorm, eps=1e-6), **new_kwargs)
    model.default_cfg = _cfg()
    return model


if __name__ == '__main__':
    vit_model = vit_base_patch2_32().cuda()
    x = torch.randn(8, 3, 128, 128).cuda()
    t = torch.randint(0, 1000, (8,)).cuda()
    c = torch.randint(0, 1000, (8,)).cuda()
    y = vit_model(x, t, c)

    print(f'x.shape {x.shape}, y.shape {y.shape}')