import torch
import pickle
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import seaborn as sns  # Import seaborn for kernel density plots


# define a custom batch processor
def batch_iter(X, batch_size, shuffle=False):
    """
    X: feature tensor (shape: num_instances x num_features)
    """
    if shuffle:
        idxs = torch.randperm(X.shape[0])
    else:
        idxs = torch.arange(X.shape[0])
    if X.is_cuda:
        idxs = idxs.cuda()
    for batch_idxs in idxs.split(batch_size):
        yield X[batch_idxs]


base_path = ''  # Define the base path for file operations.
folder = 'natl_2003_776094532_new'  # Define the folder where files will be stored.
path = os.path.join(base_path, folder, '')  # Combines the base path and folder into a complete path.
dataset_name = 'cleaned_natl2003_bin'  # Define the name of the dataset.
model_name = 'seed_776094532'

# Identify whether the system has a GPU, if yes it sets the device to "cuda:0" else "cpu"
device = "cpu" if not (torch.cuda.is_available()) else "cuda:0"

path_save = os.path.join(path, model_name)

# Load the previously saved PyTorch model from the disk
model = torch.load(path_save + '/_best_model.pt', map_location=device)

# Load original dataset
with open(path + dataset_name + '.pkl', 'rb') as f:
    data = pickle.load(f)

# Extract the DataFrame
df = data['df']
variable_list = df.columns.tolist()
# Convert DataFrame to a tensor.
df_tensor = torch.tensor(df.values, dtype=torch.float32).to(device)

# Move the model to the appropriate device (GPU if available or CPU)
model = model.to(device)

# Disable gradient computation, which is not needed and can save memory.
with torch.no_grad():
    # Set the model to evaluation mode. This deactivates dropout and batch normalization, if used.
    model.eval()

    # Perform model inference on the batched data
    z_df = []
    for cur_x in batch_iter(df_tensor, batch_size=4096, shuffle=False):
        # Pass the current batch of data through the model.
        z, _ = model(cur_x)
        z_df.append(z.cpu())  # Append the output tensor to the list, moving it to CPU

# Concatenate all batches into a single DataFrame
result_df = pd.DataFrame(torch.cat(z_df).numpy(), columns=variable_list)

# Save the DataFrame to a CSV file
result_df.to_csv(os.path.join(path, 'inverted_data.csv'), index=False)

