from torch.utils.tensorboard import SummaryWriter
import pandas as pd
import os

class TensorboardWriter(SummaryWriter):
    def __init__(self, df_save_interval:int=10, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.df = pd.DataFrame()
        self.df_update = 0
        self.df_save_interval = getattr(self, 'df_save_interval', df_save_interval)
    
    def update_df(self, data:dict):
        self.df = pd.concat([self.df, pd.DataFrame([data])], ignore_index=True)
        self.df_update += 1
        if self.df_update % self.df_save_interval == 0:
            self.save_df()
    
    def save_df(self):
        save_path = os.path.join(self.log_dir, 'tensorboard.csv')
        self.df.to_csv(save_path, index=False)