## HubRouter:
This is an implementation of the NeurIPS 2023 paper "HubRouter: Learning Global Routing via Hub Generation and Pin-hub Connection" (*HubRouter*). This approach is a global routing solver that includes a two-phase learning framework.

## Directory:
```
HubRouter/
|-- LICENSE
|-- README.md
|-- inference.py
|-- main.py
|-- main_REST.py
|-- checkpoint
|-- configs
|-- dataset
|-- logs
|-- result
|-- REST_tool
|   |-- REST_utils.py
|   |-- algorithms
|   |-- models
|   |-- test_set
|-- data
|   |-- RSMT.py
|   |-- __init__.py
|   |-- aligned.py
|   |-- condition.py
|   `-- data_interface.py
|-- models
|   |-- __init__.py
|   |-- actor_critic.py
|   |-- autoencoder.py
|   |-- ddim.py
|   |-- gan.py
|   |-- model_interface.py
|   |-- template_models.py
|   `-- vae.py
|-- preprocess
|   |-- generateInput.py
|   |-- grid.py
|   |-- test.py
|   `-- utils.py
`-- utils
    |-- __init__.py
    |-- base_utils.py
    |-- data_utils.py
    |-- diffusion_utils.py
    |-- gan_utils.py
    `-- inference_utils.py
```

## Steps
+ Add ISPD-07 benchmarks into ``./preprocess/benchmark/``;
+ Move the outputs of [NCTU](https://people.cs.nctu.edu.tw/~whliu/NCTU-GR.htm) of benchmarks into ``./preprocess/output/``;
+ Generate training dataset into ``./dataset/`` by
```
cd preprocess
python ./preprocess/generateInput.py
cd ..
```

+ Train the model in the hub-generation phase:
```
python main.py
```

+ Train the model in the pin-hub-connection phase:
```
python main_REST.py
```

+ Inference on ISPD-07 using GAN:
```
python inference.py --model GAN --case ISPD
```

## Requirements:
+ pip install -r requirements.txt

## Citations:
If you find our paper/code useful in your research, please cite
```
@inproceedings{du2023hubrouter,
  title = {HubRouter: Learning Global Routing via Hub Generation and Pin-hub Connection},
  author = {Du, Xingbo and Wang, Chonghua and Zhong, Ruizhe and Yan, Junchi},
  booktitle = {Advances in Neural Information Processing Systems},
  year = {2023}
}
```
