# README
- This is the code implementation for paper "FlexPlanner: Flexible 3D Floorplanning via Deep Reinforcement Learning in Hybrid Action Space with Multi-Modality Representation"

## Requirements
```bash
conda create -n flexplanner python=3.10
conda activate flexplanner
# CUDA 11.6
conda install pytorch==1.12.1 torchvision==0.13.1 cudatoolkit=11.6 -c pytorch -c conda-forge
# CUDA 11.3
conda install pytorch==1.12.1 torchvision==0.13.1 cudatoolkit=11.3 -c pytorch
# tensorboard and tensorboard-data-server
conda install tensorboard==2.17.0
pip install -r requirements.txt
```

## Run
```bash
bash script/alignment.sh
```
- For more arguments, please refer to `arguments/args.py`