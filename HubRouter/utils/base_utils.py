import os
import importlib
import lightning.pytorch.callbacks as plc

def instantiate_from_config(config):
    return get_obj_from_str(config["target"])(**config.get("params", dict()))

def get_obj_from_str(string, reload=False):
    module, cls = string.rsplit(".", 1)
    if reload:
        module_imp = importlib.import_module(module)
        importlib.reload(module_imp)
    return getattr(importlib.import_module(module, package=None), cls)

def examine_dir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.mkdir(dir)
        print(f'make directory: {dir}')

def load_callbacks(monitor='val_loss'):
    callbacks = []
    callbacks.append(plc.EarlyStopping(
        monitor=monitor,
        mode='min',
        patience=20
    ))

    callbacks.append(plc.ModelCheckpoint(
        monitor=monitor,
        filename='best-{epoch:02d}-{val_loss:.6f}',
        save_top_k=1,
        mode='min',
        save_last=True
    ))

    callbacks.append(plc.LearningRateMonitor(
            logging_interval='epoch'))
    return callbacks
    