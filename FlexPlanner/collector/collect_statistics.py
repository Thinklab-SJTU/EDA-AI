import tianshou
from tianshou.env import DummyVectorEnv
from tianshou.policy import BasePolicy
from tianshou.data import VectorReplayBuffer, Collector

def get_statistics(envs:DummyVectorEnv, buffer:VectorReplayBuffer, policy:BasePolicy, n_episode:int, statistics_method:int) -> dict:
    """return the mean of some metrics"""
    envs.reset()
    buffer.reset()
    collector = Collector(policy, envs, buffer)
    collector.reset()

    result = collector.collect(n_episode=n_episode)

    data, indices = buffer.sample(0)
    if statistics_method == 0:
        hpwl_mean = data.info.hpwl.mean()
    elif statistics_method == 1:
        hpwl_mean = data.info.hpwl[data.terminated].mean()
    elif statistics_method == 2:
        hpwl_mean = data.info.hpwl[data.terminated].max()
    else:
        raise NotImplementedError("statistics_method = {} is not implemented".format(statistics_method))
    
    statistics = {
        "hpwl": hpwl_mean
    }
    return statistics
