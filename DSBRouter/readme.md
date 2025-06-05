# DSBRouter: End-to-end Global Routing via Diffusion Schr¨odinger Bridge

Implementation of the paper [DSBRouter: End-to-end Global Routing via Diffusion Schr¨odinger Bridge](https://arxiv.org/abs/2403.14623).

## Training
1. Use [nthurouter](https://www.cs.nthu.edu.tw/~tcwang/nthuroute/) to produce the routing result of training datasets. Then follow the preprocess steps in [hubrouter](https://github.com/Thinklab-SJTU/EDA-AI) to produce the images of initial pins and routes.
2. Clone the repo & Initialize the env. For the torch, plz choose the right [version](https://pytorch.org/get-started/previous-versions/) according to your installed [cuda](https://developer.nvidia.com/cuda-toolkit-archive).
   
   ```bash
   conda create -n DSBRouter python=3.10
   conda activate DSBRouter
   # Get torch installed
   pip3 install torch torchvision torchaudio
   pip install -r requirements.txt
   git clone https://github.com/Thinklab-SJTU/EDA-AI.git
   cd DSBRouter
   ```

3. Prepare dataset

   Move the well-prepared pin imgaes into the folder `./DSBRouter/dataset/EDA/prior` and routes images into the folder `./DSBRouter/dataset/EDA/data`

4. Start to train

   ```
   torchrun --standalone --nproc_per_node=8 train.py \
   --method dsb --noiser flow --network uvit-b --batch_size 128 \
   --prior eda_prior --dataset eda_data --val_prior eda_prior --val_data eda_data \
   --lr 1e-5 --repeat_per_epoch 256 --use_amp --training_timesteps 50 --inference_timesteps 50 \
   --simplify --reparam term --gamma_type linear_1e-3_1e-2 --exp_name EDA_transfer
   ```


## Inference
1. Put the tested [benchmark](https://www.google.com.hk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjqgZ_7zcGNAxX5f_UHHegTEpkQFnoECAwQAQ&url=http%3A%2F%2Fwww.ispd.cc%2Fcontests%2F08%2Fispd08rc.html&usg=AOvVaw0EJgZ1S3r5OUQf0rqRs8ee&opi=89978449) into the `./DSBRouter/preprocess/benchmark` folder.
   
   ```bash
   python inference.py --network adm --prior afhq-dog-512 \
   --dataset afhq-cat-512 --simplify --reparam term \
   --gamma_type linear_1e-3_1e-2 --exp_name trdsb-afhq512 \
   --ckpt ./ckpt/afhq512.pth --num_sample 128 \
   --batch_size 16
   ```

2. Prepare dataset

   Move the well-prepared pin imgaes into the folder `./DSBRouter/dataset/EDA/prior` and routes images into the folder `./DSBRouter/dataset/EDA/data`

3. Start to train

   ```
   torchrun --standalone --nproc_per_node=8 train.py \
   --method dsb --noiser flow --network adm --batch_size 128 \
   --prior eda_prior --dataset eda_data --val_prior eda_prior --val_data eda_data \
   --lr 1e-5 --repeat_per_epoch 256 --use_amp --training_timesteps 50 --inference_timesteps 50 \
   --simplify --reparam term --gamma_type linear_1e-3_1e-2 --exp_name EDA_transfer


## Inference

Here we provide some example scripts for sampling from pre-trained models.

**AFHQ 512**

```bash
python inference.py --network adm --prior afhq-dog-512 \
   --dataset afhq-cat-512 --simplify --reparam term \
   --gamma_type linear_1e-3_1e-2 --exp_name trdsb-afhq512 \
   --ckpt ./ckpt/afhq512.pth --num_sample 128 \
   --batch_size 16
```

`--prior` sets the prior distribution ($p_{\text{prior}}$); `--dataset` is the data distribution ($p_{\text{data}}$); `--simplify` is a flag to use *Simplified DSB*; `--reparam` chooses the way for reparameterization, `term`
 means *Terminal Reparameterization*, `flow` means *Flow Reparameterization*, default is `None`; `--gamma_type` controls the way to add noise to construct $p_{\text{ref}}$; `--ckpt` points to the path of pre-trained model.

Or you could run `python inference.py -h` to see the full argument list.

**AFHQ 256**

```bash
python inference.py --network adm --prior afhq-dog-256 \
   --dataset afhq-cat-256 --simplify --reparam term \
   --gamma_type linear_1e-3_1e-2 --exp_name trdsb-afhq256 \
   --ckpt ./ckpt/afhq256.pth
```

**CelebA 64**

```bash
python inference.py --network uvit-b --prior pixel-standard \
   --dataset celeba-64 --simplify --reparam term \
   --gamma_type linear_1e-5_1e-4 --exp_name trdsb-celeba \
   --ckpt ./ckpt/celeba.pth
```

**2D experiments**

```bash
python inference_2d.py --prior dsb-pinwheel --dataset checkerboard:8 \
   --exp2d --simplify --gamma_type linear_1e-4_1e-3 \
   --exp_name sdsb-pinwheel-checkerboard8 --ckpt ./ckpt/sdsb-pinwheel-checkerboard8.pth
```

## Training

**2D experiments**

```bash
# Simplified DSB
torchrun --standalone train.py --exp2d --method dsb --prior dsb-pinwheel --dataset checkerboard:8 --training_timesteps 16 --inference_timesteps 16 --gamma_type linear_1e-4_1e-3 --repeat_per_epoch 8 --epochs 41 --exp_name sdsb-pinwheel-checkerboard --noiser flow --simplify
```

**AFHQ512**

```bash
torchrun --standalone --nproc_per_node=8 train.py --method dsb --noiser flow --network adm --batch_size 192 --prior afhq-dog-512 --dataset afhq-cat-512 --val_prior afhq-dog-512 --val_data afhq-cat-512 --lr 1e-5 --repeat_per_epoch 256 --use_amp --training_timesteps 100 --inference_timesteps 100 --simplify --reparam term --gamma_type linear_1e-3_1e-2 --exp_name trdsb-afhq512 --backward_ckpt ./ckpt/afhq512_fm_dog2cat.pth --forward_ckpt ./ckpt/afhq512_fm_cat2dog.pth --skip_epochs 1
```

For more training settings, please refer to [`training_command.md`](./training_command.md).
