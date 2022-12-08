import torch
from os.path import dirname, abspath

x = torch.load(dirname(dirname(abspath(__file__))) + '/data/dist.pt')
print(x)