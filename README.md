# DeepPlace

An end-to-end learning approach DeepPlace for placement problem with two stages. The deep reinforcement learning (DRL) agent places the macros sequentially, followed by a gradient-based optimization placer to arrange millions of standard cells. We use [PPO](https://github.com/ikostrikov/pytorch-a2c-ppo-acktr-gail) for all the experiments implemented with Pytorch, and the GPU version of
[DREAMPlace](https://github.com/limbo018/DREAMPlace) is adopted as gradient based optimization placer for arranging standard cells.

## Citations

If you find our paper/code useful in your research, please citing
```
@article{cheng2021joint,
  title={On Joint Learning for Solving Placement and Routing in Chip Design},
  author={Cheng, Ruoyu and Yan, Junchi},
  journal={Advances in Neural Information Processing Systems},
  volume={34},
  pages={16508--16519},
  year={2021}
}
```

```
@inproceedings{chengpolicy,
  title={The Policy-gradient Placement and Generative Routing Neural Networks for Chip Design},
  author={Cheng, Ruoyu and Lyu, Xianglong and Li, Yang and Ye, Junjie and Jianye, HAO and Yan, Junchi},
  booktitle={Advances in Neural Information Processing Systems}
}
```
