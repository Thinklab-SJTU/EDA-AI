# PRNet

A neural mixed-size placement and routing pipeline **PRNet**. The placement part is fulfilled by a policy gradient based RL method for macros by considering their sizes. The routing is achieved by one-shot generation of the whole path, with our devised net order learning module to dynamically adjust the routing order. We use [PPO](https://github.com/ikostrikov/pytorch-a2c-ppo-acktr-gail) for all the experiments implemented with Pytorch, and the GPU version of
[DREAMPlace](https://github.com/limbo018/DREAMPlace) is adopted as gradient based optimization placer for arranging standard cells.

## Requirements

* [PyTorch](http://pytorch.org/)
* [OpenAI baselines](https://github.com/openai/baselines)
* [DREAMPlace](https://github.com/limbo018/DREAMPlace)

In order to install requirements, follow:

```bash
# PyTorch
pip install torch==1.8.1+cu111 torchvision==0.9.1+cu111 -f https://download.pytorch.org/whl/lts/1.8/torch_lts.html

# Baselines for Atari preprocessing
git clone https://github.com/openai/baselines.git
cd baselines
pip install -e .

# Other requirements
pip install -r requirements.txt

# DREAMplace installation
git clone --recursive https://github.com/limbo018/DREAMPlace.git
mkdir build 
cd build 
cmake .. -DCMAKE_INSTALL_PREFIX=your_install_path -DPYTHON_EXECUTABLE=$(which python)
make 
make install

#Get benchmarks
python benchmarks/ispd2005_2015.py

# DGL installation
conda install -c dglteam dgl-cuda10.2
```

## Training

### Conditional Generative Router

* Train a model from scratch:
```
python train.py --dataroot $DATAROOT$ --A_folder $A_FOLDER$ --B_folder $B_FOLDER$ --model $MODEL$ --netG $NET_G$ --netD $NET_D$ --batch_size $BATCH_SIZE$ --gpu_ids $IDS$ --n_epochs $N_EPOCH$
```

We provide a couple of different models which contains purely generative models and conditional generative adversarial models, since we intend to explore the performance of diverse structures and networks.

* Continue to train:
```
python train.py --dataroot $DATAROOT$ --A_folder $A_FOLDER$ --B_folder $B_FOLDER$ --model $MODEL$ --netG $NET_G$ --netD $NET_D$ --batch_size $BATCH_SIZE$ --gpu_ids $IDS$ --n_epochs $N_EPOCH$ --epoch $EPOCH$
```

* Train the expanded model:
```
python train_expand.py --dataroot $DATAROOT$ --A_folder $A_FOLDER$ --B_folder $B_FOLDER$ --model $MODEL$ --netG $NET_G$ --netD $NET_D$ --batch_size $BATCH_SIZE$ --gpu_ids $IDS$ --n_epochs $N_EPOCH$ --load_path $PATH_OF_TRAINED_BASIC_MODEL$
```

* Test the model:
```
python test.py --eval --dataroot $DATAROOT$ --results_dir $DIR$ --A_folder $A_FOLDER$ --B_folder $B_FOLDER$ --model $MODEL$ --netG $NET_G$ --batch_size $BATCH_SIZE$ --gpu_ids $IDS$ --epoch $EPOCH$
```
The test results will be saved to an directory **DIR**, './test_results' by default.

### Overall Placement and Routing

```bash
python main.py --task "route" --algo ppo --use-gae --lr 2.5e-4 --clip-param 0.1 --value-loss-coef 0.5 --num-processes 1 --num-steps 2162 --num-mini-batch 1 --log-interval 1 --use-linear-lr-decay --entropy-coef 0.01 --log-name 3-01 --net-num 1452 --cell-num 710 --base1 28500 --base2 70000 --rate1 0.03 --rate2 0.008

```

## Results

![circuit adaptec4](imgs/adaptec4.png)

## Citations

If you find our paper/code useful in your research, please citing
```
@inproceedings{chengpolicy,
  title={The Policy-gradient Placement and Generative Routing Neural Networks for Chip Design},
  author={Cheng, Ruoyu and Lyu, Xianglong and Li, Yang and Ye, Junjie and Jianye, HAO and Yan, Junchi},
  booktitle={Advances in Neural Information Processing Systems}
}
```