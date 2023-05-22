#!/usr/bin/env python
# coding: utf-8
""" Python module for plotting connectivity results

Run as:

    python Visualize/plot_connect.py
"""

import numpy as np
import matplotlib.pyplot as plt

def get_sep_lines(labels):
    """ Get the separation lines for heatmap plot based on the network labels

    Parameters
    ----------
    labels : Pandas DataFrame
        DataFrame containing the ROI and Network name for each atlas region.

    Returns
    -------
    Pandas Series
        Series with the name of the network (as index) and the number of regions within the network.
    """

    sep_lines = labels.groupby(by = "Network", sort = False).count()

    return sep_lines.iloc[:, 0]

def plot_grid(size_between_lines, n_reg=156):
    """ Plot a grid over a heatmap

    Parameters
    ----------
    size_between_lines : np.array or list
        Indices of the space between the lines to draw
    n_reg : int, optional
        Number of regions/pixels in the heatmap, by default 156
    """

    line_idx = np.cumsum(size_between_lines)

    if line_idx.max() > n_reg:
        print("Warning, maximum line id is larger than the number of regions",
                f"({line_idx.max()} > {n_reg}) - Max values adjusted.")
        line_idx[np.where(line_idx > n_reg)[0]] = n_reg

    for line in line_idx:
        plt.plot(np.array([0, n_reg]) - .5, np.array([line]*2) - .5, c = 'silver', linewidth = 3)
        plt.plot(np.array([line]*2) - .5, np.array([0, n_reg]) - .5, c = 'silver', linewidth = 3)

    return

def heatmap_fc(data, title = '', labels = None, rot = 55, axes = None, save_loc = None, fontsize_overwrite = None, sym_cbar = True, figsize = (45, 35), **kwargs):

    PXL_SPACE = 200

    if axes is None:
        _, axes = plt.subplots(figsize = figsize)

    if labels is None:
        labels = [np.arange(data.shape[1]), np.arange(data.shape[0])]
    elif type(labels) != list:
        labels = [labels]*2

    plt.title(title, fontsize = 30)

    # Plotting data
    if sym_cbar:
        maxval = np.abs(data).max()*1.1
        minval = -maxval
    else:
        maxval = data.max()
        minval = data.min()

    heatmap = axes.imshow(data, vmin = minval, vmax = maxval, **kwargs)

    # Setting x and y ticks
    axes.set_xticks(np.arange(data.shape[1]))
    axes.set_xticks(np.arange(data.shape[1] + 1) - .5, minor = True)
    axes.set_yticks(np.arange(data.shape[0]))
    axes.set_yticks(np.arange(data.shape[0] + 1) - .5, minor = True)
    
    if fontsize_overwrite:
        adapt_fontsize = fontsize_overwrite
    else:
        adapt_fontsize = min(axes.get_figure().get_size_inches())*axes.get_figure().dpi - PXL_SPACE
        adapt_fontsize = adapt_fontsize/(np.mean([len(labels[0]), len(labels[1])]))

    axes.set_xticklabels(labels[0], fontsize = adapt_fontsize)
    axes.set_yticklabels(labels[1], fontsize = adapt_fontsize)

    # Rotating tick labels and setting alignments
    plt.setp(axes.get_xticklabels(), rotation = rot, ha = "right", rotation_mode = "anchor")

    #plt.grid(alpha = .5)

    cbar = plt.colorbar(heatmap, aspect=50, pad = .01)
    cbar.ax.tick_params(labelsize = 60)

    plt.tight_layout()

    if save_loc is not None:
        plt.savefig(save_loc)

    return heatmap