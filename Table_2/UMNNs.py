import numpy as np
np.set_printoptions(precision=3, suppress=None)
import os
os.environ["CUDA_VISIBLE_DEVICES"] = '0'
import collections.abc
collections.Iterable = collections.abc.Iterable
from medflow import train_med, sim_med
import time

start_time = time.time()

base_path = '/project/wodtke/cGNF_python_code'  # Define the base path for file operations.
folder = 'natl_2003'  # Define the folder where files will be stored.
path = os.path.join(base_path, folder, '')  # Combines the base path and folder into a complete path.
dataset_name = 'cleaned_natl2003_bin'  # Define the name of the dataset.

if not (os.path.isdir(path)):  # checks if a directory with the name 'path' exists.
    os.makedirs(path)  # if not, creates a new directory with this name. This is where the logs and model weights will be saved.

## MODEL TRAINING
train_med(path=path, dataset_name=dataset_name, treatment='a', confounder=["c1", "c2", "c3", "c4"], mediator=["l", "m"], outcome='y',
           test_size=0.2, cat_var=["a", "l", "m", "y", "c1", "c2", "c3", "c4"], sens_corr=None, seed_split=1,
           model_name=path + 'seed_1',
           trn_batch_size=128, val_batch_size=4096, learning_rate=1e-4, seed=1,
           nb_epoch=50000, nb_estop=50, val_freq=1,
           emb_net=[100, 90, 80, 70, 60],
           int_net=[60, 50, 40, 30, 20])

## EFFECT ESTIMATION
# Path-specific effects
sim_med(path=path, dataset_name=dataset_name, model_name=path + 'seed_1', n_mce_samples=100000, seed=1, inv_datafile_name='1_path_100k',
        cat_list=[0, 1], moderator=None)

# Interventional effects
sim_med(path=path, dataset_name=dataset_name, model_name=path + 'seed_1', n_mce_samples=100000, intv_med=["m=intv"], seed=1, inv_datafile_name='1_intv_100k',
        cat_list=[0, 1], moderator=None)

end_time = time.time()
time_taken = end_time - start_time

print(time_taken)
