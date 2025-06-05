# DSBRouter: End-to-end Global Routing via Diffusion Schr¨odinger Bridge

Implementation of the paper [DSBRouter: End-to-end Global Routing via Diffusion Schr¨odinger Bridge](https://arxiv.org/abs/2403.14623). **The overall training logic code comes from https://github.com/chrisway613/SDSB**
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
   # --benchmark can be chosen in ['ibm-ada','abl-NN','abl-RST','IPSD-s4','IPSD-s','IPSD-b4','IPSD-b'] to reproduce the results in the paper.
   python inference.py --network uvit-b --benchmark ibm-ada \
   --simplify --reparam term --gamma_type linear_1e-3_1e-2 \
   --exp_name EDA_infer --ckpt $ckpt path$ --num_sample 128 \
   --batch_size 16
   ```
