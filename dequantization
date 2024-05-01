import pandas as pd
import numpy as np
import os

# Setting seeds
np.random.seed(42)

base_path = ''  # Define the base path for file operations.
folder = 'natl_2003_549243471'  # Define the folder where files will be stored.
path = os.path.join(base_path, folder, '')

# Read in the data
df = pd.read_csv(path + 'cleaned_natl2003_bin.csv')

# Create a new DataFrame to store the transformed data
transformed_df = pd.DataFrame(index=df.index)

# Iterate over each column in the DataFrame
for column in df.columns:
    unique_values = df[column].unique()  # Get unique values in the column
    column_data = np.empty(df.shape[0])  # Empty array to store transformed data

    # Process each unique value
    for value in unique_values:
        # Select indices for rows with the current category
        indices = df[column] == value
        # Draw random samples from a normal distribution centered at the category value with a standard deviation of 1/6
        column_data[indices] = np.random.normal(loc=value, scale=1/8, size=np.sum(indices))

    # Assign the transformed data to the new DataFrame
    transformed_df[column] = column_data

# Save the transformed DataFrame to a new CSV file
transformed_df.to_csv(path + "dequantized_data.csv", index=False)

