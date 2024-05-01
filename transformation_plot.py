import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import seaborn as sns
import os

plt.rcParams["font.family"] = "Times New Roman"

# Path and data setup
base_path = ''  # Define the base path for file operations.
folder = 'natl_2003_549243471'  # Define the folder where files will be stored.
path = os.path.join(base_path, folder, '')  # Combines the base path and folder into a complete path.

# Load dataframes
original_df = pd.read_csv(path + 'cleaned_natl2003_bin.csv')
dequantized_df = pd.read_csv(path + 'dequantized_data.csv')
inverted_df = pd.read_csv(path + 'inverted_data.csv')

# Define columns to plot
columns_to_plot = ['l', 'm', 'y']

# Setup plot grid
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))  # Adjust figsize as needed
fig.subplots_adjust(hspace=0.4, wspace=0.4)  # Adjust spacing between plots

data_frames = [original_df, dequantized_df, inverted_df]

titles = ['Original Data', 'Dequantized Data', 'Inverted Data']
column_header = ['$L$', '$M$', '$Y$']

# Set column headers
for i, col in enumerate(column_header):
    axs[0, i].set_title(col, fontsize=16, weight='bold')

# Set row titles
for i, row in enumerate(titles):
    axs[i, 0].text(-0.3, 0.5, row, rotation=90, verticalalignment='center', horizontalalignment='right', transform=axs[i, 0].transAxes, fontsize=12, weight='bold')


# Iterate over each type of data (row) and each feature (column)
for i, df in enumerate(data_frames):
    for j, column in enumerate(columns_to_plot):
        ax = axs[i, j]
        data = df[column.strip('$')]

        if i == 2:  # Only for inverted data
            count, bins, ignored = ax.hist(data, bins=200, density=True, alpha=0.6, color='gray')
            x = np.linspace(min(bins), max(bins), 1000)
            p = norm.pdf(x, 0, 1)
            ax.plot(x, p, 'k', linewidth=2)

            ax.set_ylabel('Density')

        elif i == 1:  # Dequantized data

            # count, bins, ignored = ax.hist(data, bins=200, density=True, alpha=0.6, color='gray')

            counts, bin_edges = np.histogram(data, bins=200, density=True)
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            adjusted_counts = counts * np.diff(bin_edges) * 20
            ax.bar(bin_centers, adjusted_counts, width=np.diff(bin_edges), alpha=0.6, color='gray')

            ax.set_ylabel('Density')
            # ax.clear()
            # adjusted_count = count * np.diff(bins)
            # bin_centers = 0.5 * (bins[:-1] + bins[1:])
            # ax.bar(bin_centers, adjusted_count, width=np.diff(bins), alpha=0.6, color='gray')

        else:  # Original data

            # Define bins such that 0 and 1 are centered
            bin_edges = np.array([-0.5, 0.5, 1.5])  # Defines bins with 0 and 1 as centers
            counts, _ = np.histogram(data, bins=bin_edges)
            bin_centers = [0, 1]  # Directly setting bin centers to 0 and 1
            total_data_points = len(data)
            adjusted_counts = counts / (total_data_points * np.diff(bin_edges))

            if column != 'm':

                # Plot droplines and circles at the top of each dropline
                for center, count in zip(bin_centers, adjusted_counts):
                    ax.vlines(center, 0, count, color='black', linewidth=0.5)  # Droplines from count to axis
                    ax.scatter(center, count, color='black', s=10, zorder=3)  # Solid circle at the top of the line

            else:  # for column 'm'
                adjusted_counts[1] *= 100  # Scale up the probability for '1'

                # Plot droplines for 0 and scaled line for 1
                for center, count in zip([0, 1], adjusted_counts):
                    if center == 0:
                        ax.vlines(center, 0.53, count, color='black', linewidth=0.5)  # Droplines from count to upper breakage
                        ax.vlines(center, 0.47, 0.53, color='grey', linestyles='dotted', linewidth=0.5)  # Droplines from upper breakage to lower breakage
                        ax.vlines(center, 0, 0.47, color='black', linewidth=0.5)  # Droplines from lower breakage to axis
                        ax.scatter(center, count, color='black', s=10, zorder=3)  # Solid circle at the top of the line
                    else:
                        ax.vlines(center, 0, count, color='black', linewidth=0.5)  # Droplines from count to axis
                        ax.scatter(center, count, color='black', s=10, zorder=3)  # Solid circle at the top of the line

                    # Add broken axis effect between 0.4 and 0.6
                    d = .01  # Proportion of vertical to horizontal extent of the diagonal lines
                    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
                    ax.plot((-d, +d), (0.48, 0.46), **kwargs)  # Lower diagonal on y axis
                    ax.plot((-d, +d), (0.54, 0.52), **kwargs)  # Upper diagonal on y axis

                    ax.spines['left'].set_bounds(-0.05, 0.47)
                    ax2 = ax.twinx()
                    ax2.spines['left'].set_bounds(0.53, 1)


                    # Manually setting y-axis ticks and labels
                    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
                    ax.set_yticklabels(['0', '0.0002', '0.0004', '0.6', '0.8', '1'])

            ax.set_ylabel('Probability')
            ax.set_xticks(bin_centers)
            ax.set_xticklabels(['0', '1'])

        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        # ax.clear()
            # adjusted_count = count * np.diff(bins)
            # bin_centers = 0.5 * (bins[:-1] + bins[1:])
            # ax.bar(bin_centers, adjusted_count, width=np.diff(bins), alpha=0.6, color='gray')

# plt.tight_layout()

# Save the entire plot
plt.savefig(f"{path}/combined_plots.png", dpi=300)
plt.close()  # Close the plot to free memory
# plt.show()  # Optionally display the plot
