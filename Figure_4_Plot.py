import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os

plt.rcParams["font.family"] = "Times New Roman"

# Path and data setup
base_path = ''  # Define the base path for file operations.
folder = 'natl_2003_250470397'
path = os.path.join(base_path, folder, '')

# Load dataframes
original_df = pd.read_csv(path + 'cleaned_natl2003_bin.csv')
dequantized_df = pd.read_csv(path + 'dequantized_data.csv')
inverted_df = pd.read_csv(path + 'inverted_data.csv')

# Define columns to plot
columns_to_plot = ['l', 'm', 'y']

# Setup plot grid
fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(15, 10))
fig.subplots_adjust(hspace=0.4, wspace=0.4)

data_frames = [original_df, dequantized_df, inverted_df]
titles = ['Original Data', 'Dequantized Data', 'Transformed Data']
column_header = ['$L$', '$X$', '$Y$']

# Set column headers
for i, col in enumerate(column_header):
    axs[0, i].set_title(col, fontsize=16, weight='bold')

# Set row titles
for i, row in enumerate(titles):
    axs[i, 0].text(
        -0.3, 0.5, row,
        rotation=90,
        verticalalignment='center',
        horizontalalignment='right',
        transform=axs[i, 0].transAxes,
        fontsize=12, weight='bold'
    )

# Plotting loop
for i, df in enumerate(data_frames):
    for j, column in enumerate(columns_to_plot):
        ax = axs[i, j]
        data = df[column]

        if i == 2:
            # Transformed data: histogram + N(0,1)
            _, bins, _ = ax.hist(data, bins=200, density=True, alpha=0.6, color='gray')
            x = np.linspace(bins.min(), bins.max(), 1000)
            ax.plot(x, norm.pdf(x, 0, 1), 'k', linewidth=2)
            ax.set_ylabel('Density')

        elif i == 1:
            # Dequantized data
            if column == 'l':
                counts, edges = np.histogram(data, bins=400, density=True)
                centers = 0.5 * (edges[:-1] + edges[1:])
                heights = counts * np.diff(edges) * 20
                ax.bar(centers, heights, width=np.diff(edges), alpha=0.6, color='gray')
                ax.set_ylabel('Density')
                ax.set_xlim(-0.6, 5.6)
                ax.set_ylim(0, 1)
                # inset zoom on l=4 and l=5
                inset_l = inset_axes(ax, width="30%", height="30%", loc='upper right', borderpad=1)
                inset_l.bar(centers, heights, width=np.diff(edges), alpha=0.6, color='gray')
                inset_l.set_xlim(3.4, 5.6)
                mask_l = (centers >= 3.4) & (centers <= 5.6)
                ymax_l = heights[mask_l].max() * 1.1
                inset_l.set_ylim(0, ymax_l)
                # ticks for l zoom
                xticks_l = [3.5, 4, 4.5, 5, 5.5]
                inset_l.set_xticks(xticks_l)
                inset_l.set_xticklabels([f"{x:.1f}" for x in xticks_l], fontsize=8)
                yticks_l = np.linspace(0, ymax_l, 3)
                inset_l.set_yticks(yticks_l)
                inset_l.set_yticklabels([f"{y:.3f}" for y in yticks_l], fontsize=8)

                # arrow indicating zoom
                idx1 = np.argmin(np.abs(centers - 1))
                y1 = heights[idx1]
                ax.annotate(
                    '',
                    xy=(4.5, y1-0.02), xycoords='data',
                    xytext=(0.85, 0.8), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', color='black', lw=1)
                )

            elif column == 'm':
                counts, edges = np.histogram(data, bins=200, density=True)
                centers = 0.5 * (edges[:-1] + edges[1:])
                heights = counts * np.diff(edges) * 20
                ax.bar(centers, heights, width=np.diff(edges), alpha=0.6, color='gray')
                ax.set_ylabel('Density')
                ax.set_xlim(-0.6, 1.6)
                ax.set_ylim(0, 0.675)
                # inset zoom on m ≈ 1
                inset = inset_axes(ax, width="30%", height="30%", loc='upper right', borderpad=1)
                inset.bar(centers, heights, width=np.diff(edges), alpha=0.6, color='gray')
                # adjust x-range to include more bins
                inset.set_xlim(0.5, 1.5)
                mask = (centers >= 0.5) & (centers <= 1.5)

                # auto-scale to highest spike in window
                ymax = heights[mask].max() * 1.1
                inset.set_ylim(0, ymax)

                # set ticks dynamically
                xticks_m = [0.6, 1, 1.4]
                inset.set_xticks(xticks_m)
                inset.set_xticklabels([f"{x:.1f}" for x in xticks_m])
                yticks = np.linspace(0, ymax, 3)
                inset.set_yticks(yticks)
                inset.set_yticklabels([f"{y:.3f}" for y in yticks])

                # arrow indicating zoom
                idx1 = np.argmin(np.abs(centers - 1))
                y1 = heights[idx1]
                ax.annotate(
                    '',
                    xy=(1, y1), xycoords='data',
                    xytext=(0.77, 0.8), textcoords='axes fraction',
                    arrowprops=dict(arrowstyle='->', color='black', lw=1)
                )

            else:
                counts, edges = np.histogram(data, bins=200, density=True)
                centers = 0.5 * (edges[:-1] + edges[1:])
                heights = counts * np.diff(edges) * 20
                ax.bar(centers, heights, width=np.diff(edges), alpha=0.6, color='gray')
                ax.set_ylabel('Density')
                ax.set_xlim(-0.6, 1.6)
                ax.set_ylim(0, 0.675)

        else:
            # Original data
            if column == 'l':
                vals = sorted(data.dropna().unique())
                edges = np.arange(vals[0] - 0.5, vals[-1] + 1.5, 1)
                counts, _ = np.histogram(data, bins=edges)
                centers = vals
                props = counts / (len(data) * np.diff(edges))
                for c, p in zip(centers, props):
                    ax.vlines(c, 0, p, color='black', linewidth=0.5)
                    ax.scatter(c, p, color='black', s=10, zorder=3)
                ax.set_xlim(-0.6, 5.6)
                ax.set_ylim(0, props.max() * 1.1)
                ax.set_ylabel('Probability')

            elif column == 'm':
                edges = np.array([-0.5, 0.5, 1.5])
                counts, _ = np.histogram(data, bins=edges)
                centers = [0, 1]
                props = counts / (len(data) * np.diff(edges))
                props[1] *= 100
                for c, p in zip(centers, props):
                    if c == 0:
                        ax.vlines(c, 0, 0.47, color='black', linewidth=0.5)
                        ax.vlines(c, 0.53, p, color='black', linewidth=0.5)
                        ax.scatter(c, p, color='black', s=10, zorder=3)
                    else:
                        ax.vlines(c, 0, p, color='black', linewidth=0.5)
                        ax.scatter(c, p, color='black', s=10, zorder=3)
                d = .01
                kw = dict(transform=ax.transAxes, color='k', clip_on=False)
                ax.plot((-d, +d), (0.48, 0.46), **kw)
                ax.plot((-d, +d), (0.54, 0.52), **kw)
                ax.spines['left'].set_bounds(-0.05, 0.49)
                ax2 = ax.twinx()
                ax2.spines['left'].set_bounds(0.53, 1)
                ax.set_yticks([0, 0.2, 0.4, 0.6, 1])
                ax.set_yticklabels(['0', '0.0002', '0.0004', '0.9', '1.0'])
                ax.set_xlim(-0.6, 1.6)
                ax.set_ylim(bottom=0)
                ax.set_ylabel('Probability')

            else:
                edges = np.array([-0.5, 0.5, 1.5])
                counts, _ = np.histogram(data, bins=edges)
                centers = [0, 1]
                props = counts / (len(data) * np.diff(edges))
                for c, p in zip(centers, props):
                    ax.vlines(c, 0, p, color='black', linewidth=0.5)
                    ax.scatter(c, p, color='black', s=10, zorder=3)
                ax.set_xlim(-0.6, 1.6)
                ax.set_ylim(0, 1.05)
                ax.set_ylabel('Probability')

        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

# Save and close
plt.savefig(os.path.join(path, 'combined_plots.png'), dpi=300)
plt.close()
