import torch
from torch import nn
from torch.distributions import Normal
from einops.layers.torch import Rearrange


class VanillaNormal(Normal):
    def __init__(self, feat:dict[str, torch.Tensor]) -> None:
        super().__init__(feat['mean'], feat['std'])

    def sample(self, sample_shape=torch.Size()) -> torch.Tensor:
        sample = super().sample(sample_shape)
        return sample

    def rsample(self, sample_shape=torch.Size()) -> torch.Tensor:
        sample = super().rsample(sample_shape)
        return sample

    def __repr__(self) -> str:
        res = "{}(mean={}, std={})".format(self.__class__.__name__, self.mean, self.stddev)
        return res
    

class RatioDecider(nn.Module):
    def __init__(self, input_dim:int, ratio_range:list[float,float], ratio_area_in_dim:int, hidden_dim:int, 
                 share_with_critics:bool, input_next_block:int):
        """
        @ratio_area_in_dim: if > 0, concat area information to the output of the backbone. It does not control the input, it is just a binary flag.
        """
        super().__init__()
        self.input_next_block = input_next_block
        
        if share_with_critics: # use the input from shared encoder to decide the layer
            print("[INFO] Use shared encoder in RatioDecider")
            self.backbone = nn.Sequential(
                nn.Linear(input_dim, hidden_dim),
                nn.ReLU(),
                nn.Linear(hidden_dim, hidden_dim//2),
                nn.ReLU(),
            )
        else: # use original stacked mask to decide the layer
            self.backbone = nn.Sequential(
                nn.Conv2d(input_dim, hidden_dim, 3, 1, 1),
                nn.ReLU(),
                nn.Conv2d(hidden_dim, hidden_dim, 3, 1, 1),
                nn.ReLU(),
                nn.Conv2d(hidden_dim, hidden_dim//2, 3, 1, 1),
                nn.AdaptiveMaxPool2d((1, 1)),
                Rearrange('b c 1 1 -> b c'),
            )

            
        self.share_with_critics = share_with_critics

        self.ratio_area_in_dim = max(0, ratio_area_in_dim)
        if self.ratio_area_in_dim > 0:
            self.fc_area = nn.Identity()
            fc_area_dim = 1
        else:
            fc_area_dim = 0

        self.fc_mean = nn.Sequential(
            nn.Linear(hidden_dim//2 + fc_area_dim, 1),
            nn.Tanh(),
        )


        self.logstd = nn.Parameter(torch.zeros(1))

        self.ratio_range = ratio_range
        if self.ratio_range is not None:
            self.low, self.high = self.ratio_range


    def forward(self, output_of_shared_encoder:torch.Tensor, stacked_mask:torch.Tensor, grid_area_next:torch.Tensor) -> torch.Tensor:
        """
        Input:
            x: [B, C] tensor
            grid_area_next: [B,] tensor
        Output: should return mean and std without any linear transformation
            mean: [B,] tensor
            std: [B,] tensor
        """
        if self.share_with_critics:
            x = output_of_shared_encoder
        else:
            x = stacked_mask

        x = self.backbone(x)

        if self.ratio_area_in_dim > 0: # concat area information
            grid_area_feat = grid_area_next.unsqueeze(-1)
            grid_area_feat = self.fc_area(grid_area_feat)
            x = torch.cat([x, grid_area_feat], dim=-1)

        mean = self.fc_mean(x).squeeze(-1)
        logstd = self.logstd.expand_as(mean)
        std = logstd.exp()

        return mean, std
    
