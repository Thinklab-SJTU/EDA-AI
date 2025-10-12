# OAREST README

This is an implementation of the paper **Train on Pins and Test on Obstacles for Rectilinear Steiner Minimum Tree** (OAREST).  

![OAREST Pipeline.](https://anonymous.4open.science/api/repo/OAREST_figures-B8BE/file/OAREST_pipeline.png?v=17bee186)

*Figure 1*: **OAREST Pipeline**. a) Actor network that generates the rectilinear edge sequence (RES) step-by-step with dynamic masking strategies. b) Critic network that fits the actual length of RES. c) Dynamic masking strategies, including input masking, obstacle masking, visited masking, and activation masking, that enable the RES generation with multi-degree forward passing and obstacle avoidance. d) Visualization of the connecting process of OAREST on a toy sample with 3 pins ($p_0, p_1, p_2$) and 3 obstacles ($b_0, b_1, b_2$). 

## Prerequisites
- Please run `pip install -r requirements.txt` to achieve the environment.
- This repo is executed under `torch=2.4.0+cu118` and `pytorch-lightning=2.4.0`. Please find the suitable versions of [torch](https://pytorch.org/) and [pytorch-lightning](https://lightning.ai/docs/pytorch/stable/versioning.html#compatibility-matrix).

## Datasets
We use randomly generated test data and real-world global routing benchmarks from ICCAD19.
#### Randomly Generated Test Data
[R5-R50](https://github.com/cuhk-eda/REST/tree/main/test_set). This dataset includes varying degrees ranging from 5 to 50 in increments of 5, referred to as R5, R10, etc. Each subset contains 10k test instances.
#### ICCAD19
[ICCAD19](https://www.iccad-contest.org/2019/problems.html). We use the global routing benchmarks in Problem C.

## Usage
#### Train
```
python main.py -d OAREST_train
```

#### Inference
```
python main.py -d OAREST_test
```

## Visualizations
![RSMT.](https://anonymous.4open.science/api/repo/OAREST_figures-B8BE/file/r50_o0_gst-5.262_oarest-5.265.png?v=a2e99327)
*Figure 2*: RSMT with 50 pins.
![OARSMT.](https://anonymous.4open.science/api/repo/OAREST_figures-B8BE/file/r20_o10_gst-3.551_oarest-3.845.png?v=cc95a382)
*Figure 3*: OARSMT with 20 pins and 10 obstacles.

## Citations:
If you find our paper/code useful in your research, please cite
```
@inproceedings{du2025OAREST,
  title={Train on Pins and Test on Obstacles for Rectilinear Steiner Minimum Tree},
  author={Xingbo Du, Ruizhe Zhong, Junchi Yan},
  booktitle = {Advances in Neural Information Processing Systems},
  year = {2025}
}
```

## Directory

```
|-- README.md
|-- main.py
|-- requirements.txt
|-- REST_tool
    |-- REST_utils.py
    |-- algorithms
|-- configs
    |-- OAREST_test.yaml
    |-- OAREST_train.yaml
|-- data
    |-- __init__.py
    |-- RSMT.py
    |-- data_interface.py
|-- images
    |-- OAREST_pipeline.png
|-- models
    |-- __init__.py
    |-- actor_critic.py
    |-- model_interface.py
    |-- self_attn.py
    |-- template_nets.py
|-- utils
    |-- __init__.py
    |-- base_utils.py
    |-- obstacle_utils.py
    |-- plot_utils.py
    |-- rsmt_utils.py
```
## License

OAREST is released under the **CU-SD License** (a BSD-like open license).

It includes code derived from the [REST](https://github.com/cuhk-eda/REST) project  
Copyright (c) 2022, The Chinese University of Hong Kong  
used under the CU-SD License.

Modifications and additional work © 2025 Shanghai Jiao Tong University.  
See the [LICENSE](./LICENSE) file for full details.