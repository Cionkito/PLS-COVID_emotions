#!/usr/bin/env python
# coding: utf-8
""" Python module for plotting results from the Partial Least Square Correlation analyses

Run as:

    python Visualize/plot_PLSC.py
"""

import numpy as np
import pandas as pd
import os.path as op
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors
import h5py
import pycircos

import sys

sys.path.append(op.abspath(f'../'))
#print(op.abspath('../'))

import plot_connect as vpc
import Utility.fmri_utils as futils


def load_PLSC_image(study, path_to_study, LC_id = 0, results = "BSR", n_area = 156):
    """Load the data from a PLSC analysis conducted with mypls.

    Parameters
    ----------
    study : str
        Name of the study to load
    path_to_study : str
        Path to the study mypls directory
    LC_id : int, optional
        Index of the latent component, by default 0
    results : str, optional
        Type of results to load ("BSR", "behav_loading", "img_loading", etc), by default "BSR"
    n_area : int, optional
        Number of area in the fMRI atlas/parcellation, by default 156

    Returns
    -------
    Numpy array
        PLSC data
    """

    NORM = 1
    if "GRP" in study:
        NORM = 2
        
    path_to_file = op.join(op.realpath(path_to_study),study,
                            f'myPLS_behavior_norm{NORM:1d}-{NORM:1d}_res.mat')

    try:
        f = h5py.File(path_to_file, 'r')
    except FileNotFoundError:
        print(f"File not found with specified normalization. Trying {3-NORM:1d} instead.")
        path_to_file = op.join(op.realpath(path_to_study),study,
                                f'myPLS_behavior_norm{3-NORM:1d}-{3-NORM:1d}_res.mat')
        
        f = h5py.File(path_to_file, 'r')


    if results == "BSR":
        # Computation of Bootrstrap Ratio
        PLS_image = f['res']['V'][LC_id]/f['res']['boot_results']['Vb_std'][LC_id]
    elif results == "img_loadings":
        PLS_image = f['res']['LC_img_loadings'][LC_id]
    else:
        try:
            PLS_image = f['res'][results][LC_id]
        except:
            print(f'Unknown "myPLS" results field: {results}')
            PLS_image = f['res'][results][LC_id]

    PLS_image[np.isnan(PLS_image)] = 0

    image_array = np.zeros((n_area, n_area))

    triu_idx = np.triu_indices_from(image_array, k = 1)

    # Storing data vector in a matrix
    image_array[triu_idx[0], triu_idx[1]] = PLS_image
    image_array[triu_idx[1], triu_idx[0]] = PLS_image

    f.close()

    return image_array

def plot_imaging(study, path_to_study, LC_id = 0, results = "BSR", threshold = 0, by_network = False, roi = None,
    path_to_atlas = "D:\cioncaa\Documents\Parcellations\Custom_Yeo_WholeBrain", axes = None, **kwargs):
    """ Plot the imaging pattern of a PLSC analysis

    Parameters
    ----------
    study : str
        Name of the study
    path_to_study : str
        Path to the study directory 
    LC_id : int, optional
        Index of the latent component, by default 0
    results : str, optional
        Type of the results to analyse ("loadings", "BSR", etc), by default "BSR"
    threshold : int, optional
        Threshold to the image, by default 0
    by_network : bool, optional
        Condition to sort the regions by network and draw lines between networks, by default False
    path_to_atlas : str, optional
        Path to the atlas directory, by default "D:\cioncaa\Documents\Parcellations\Custom_Yeo_WholeBrain"
    axes : matplotlib.axes, optional
        Axis to draw the figure, by default None
    """

    if axes is None:
        _, axes = plt.subplots(figsize = (45, 35))

    if roi is None:
        roi = np.arange(156)
    
    n_area = len(roi)

    plsc_data = load_PLSC_image(study, path_to_study, LC_id=LC_id, results=results, n_area=n_area)

    atlas_labels = futils.get_atlas_labels(path_to_atlas, filters = "all")

    atlas_labels = atlas_labels.iloc[roi]

    if by_network:
        atlas_labels, sorting_array = futils.sort_atlas_labels(atlas_labels)
        plsc_data = plsc_data[sorting_array][:, sorting_array]

    to_plot = plsc_data

    if threshold:
        plsc_data[np.abs(plsc_data) < threshold] = 0
        survivors = np.unique(np.where(plsc_data != 0)[0])
        to_plot = plsc_data[survivors][:, survivors]
        atlas_labels = atlas_labels.iloc[survivors]

    vpc.heatmap_fc(to_plot, axes=axes, labels=atlas_labels.roi_name, **kwargs)

    if by_network:
        sep_lines = vpc.get_sep_lines(atlas_labels)
        vpc.plot_grid(sep_lines, to_plot.shape[0])

    return to_plot, atlas_labels#, sep_lines

def get_circus_data(data, labels, net_visualize = False):

    circus_data = pd.DataFrame(index = np.arange(len(np.triu_indices_from(data, k = 1)[0])),
                                columns = ["net_from", "id_from_start", "net_to", "id_to_start", "strength",
                                            "reg_from", "reg_to", "score"])

    networks = vpc.get_sep_lines(labels)
    offset = np.cumsum(networks).shift(1, fill_value = 0)

    line_idx = 0
    maxval = np.abs(data).max()
    for from_idx in range(data.shape[0]):
    #for from_idx in [0, 1, 2]:
        for to_idx in range(from_idx+1, data.shape[0], 1):
            if np.abs(data[from_idx, to_idx]) > 0:

                net_from = labels.iloc[from_idx, -2]
                circus_data.loc[line_idx, "net_from"] = net_from
                circus_data.loc[line_idx, "id_from_start"] = from_idx - offset[net_from]
                
                net_to = labels.iloc[to_idx, -2]
                circus_data.loc[line_idx, "net_to"] = net_to
                circus_data.loc[line_idx, "id_to_start"] = to_idx - offset[net_to]

                strength_percent = 50*data[from_idx, to_idx]/maxval+50
                circus_data.loc[line_idx, "strength"] = strength_percent

                circus_data.loc[line_idx, "reg_from"] = labels.iloc[from_idx, 0]
                circus_data.loc[line_idx, "reg_to"] = labels.iloc[to_idx, 0]
                circus_data.loc[line_idx, "score"] = data[from_idx, to_idx]

                line_idx += 1

    if net_visualize:
        net_from, net_to = [], []

        for i, row in networks.iteritems():
            net_to += networks[i:].index.tolist()
            net_from += [i]*len(networks[i:])

        id_from, id_to = [], []
        for i in range(19):
            id_from.append(0)
            id_from += np.arange(18, i, -1).tolist()

            id_to.append(19)
            id_to += (np.arange(18-i)+1).tolist()

        circus_by_net = pd.DataFrame(np.array([net_from, net_to]).T,
                                columns = ["net_from", "net_to"])
        circus_grouped = circus_data.groupby(by = ["net_from", "net_to"], sort = False).strength.aggregate(['count', 'mean']).reset_index()#.mean().reset_index()

        circus_by_net = circus_by_net.merge(circus_grouped, left_on=['net_from', 'net_to'],
                    right_on=['net_from', 'net_to'], how = "left")

        circus_by_net.loc[:, "id_from_start"] = id_from
        circus_by_net.loc[:, "id_to_start"] = id_to

        circus_data = circus_by_net.rename(columns={"count":"alpha",
                                        "mean":"strength"}).fillna(0)

        circus_data.alpha = circus_data.alpha.apply(lambda x: x/circus_data.alpha.max())

        #circus_data = circus_by_net.fillna(0)

    return circus_data.dropna(axis = 0)

def get_stronger_fc(data, labels, n_val = 15, use_abs = False, **kwargs):

    circus_data = get_circus_data(data, labels)

    if use_abs:
        sorted_circus = circus_data.sort_values(by='score', key = np.abs, **kwargs)
    else:
        sorted_circus = circus_data.sort_values(by='score', **kwargs)

    sorted_circus.loc[:, 'reg_from'] = sorted_circus.reg_from.apply(futils.get_roi_name)
    sorted_circus.loc[:, 'reg_to'] = sorted_circus.reg_to.apply(futils.get_roi_name)

    return sorted_circus.iloc[:n_val, [0, 5, 2, 6, 7]]

def plot_circus(data, labels, fig = None, figsize = (15, 15), net_visualize = False, label_size=30, **kwargs):
    Garc = pycircos.Garc
    Gcircle = pycircos.Gcircle

    networks = vpc.get_sep_lines(labels)

    circus_data = get_circus_data(data, labels, net_visualize=net_visualize)

    if fig is None:
        circle = Gcircle(figsize=figsize)
    else:
        circle = Gcircle(fig=fig, figsize=fig.get_size_inches())
    arc_cm = cm.get_cmap('Spectral', 21)

    MULTLINK = 1
    STRENGTH_PERCENTILE = 85
    ARC_RADIUS = 800#900

    if net_visualize:
        for i, (net, size) in enumerate(zip(networks.index, networks)):
            arc = Garc(arc_id=net, size=len(networks)+2, interspace=2,
                       raxis_range=(ARC_RADIUS, ARC_RADIUS+50), labelposition=80,
                       label_visible=True, labelsize=label_size,
                       facecolor=colors.to_rgba(arc_cm(i+1)))
            circle.add_garc(arc)
    else:
        MULTLINK = 10
        for i, (net, size) in enumerate(zip(networks.index, networks)):
            arc = Garc(arc_id=net, size=size*MULTLINK+1, interspace=2,
                       raxis_range=(ARC_RADIUS, ARC_RADIUS+50), labelposition=45,
                       label_visible=True, labelsize=label_size,
                       facecolor=colors.to_rgba(arc_cm(i+1)))
            circle.add_garc(arc)

    circle.set_garcs()

    links_cm = cm.get_cmap('coolwarm', 100)

    if net_visualize:
        strength_thresh = np.percentile(np.abs(circus_data.strength), STRENGTH_PERCENTILE, interpolation = 'nearest')
        #circus_data['alpha'] = circus_data.strength.apply(lambda x: 0.5*np.abs(x)/strength_thresh)
    else:
        circus_data['alpha'] = circus_data.strength.apply(lambda x: 0.1*np.abs(x-50)/50)

    for row_idx, row in circus_data.iterrows():
        if np.abs(row.strength) > 0:
            source = (row.net_from, row.id_from_start*MULTLINK, (row.id_from_start+1)*MULTLINK, ARC_RADIUS+50)
            destin = (row.net_to, row.id_to_start*MULTLINK, (row.id_to_start+1)*MULTLINK, ARC_RADIUS+50)
            circle.chord_plot(source, destin, colors.to_rgba(links_cm(int(row.strength), alpha = row.alpha)))

    return

def plot_bar(data, labels, scatter = None, groupings = ['Mild', 'Moderate', 'Severe'], colors = ['tab:green', 'tab:orange', 'tab:red'],
            rot = 0, HIGHLIGHT = True, fig = None, ax = None, figsize = (35, 20), ylabels = 'Correlations', **kwargs):
    nlab = len(labels)

    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize = figsize)

    offsets = np.linspace(-0.4, 0.4, len(groupings)+2)[1:-1]
    if offsets[0]:
        width = offsets[0]
    else:
        width = .5

    for i, (x_off, col, grp) in enumerate(zip(offsets, colors, groupings)):
    #for i, (x_off, col, grp) in enumerate(zip([-.25, 0, .25], colors, groupings)):
    #for i, (x_off, col) in enumerate(zip([-.2, .2], ['tab:red', 'tab:blue'])):
        x_loc = np.arange(nlab)+x_off
        bar_height = data.iloc[nlab*i:nlab*(i+1), 0]
        plt.bar(x_loc, bar_height, width=width,
                color = col, alpha = .7, edgecolor = 'k', linewidth = 4, label = grp, **kwargs)

        if scatter is not None:
            for j in range(nlab):
                plt.scatter(np.random.normal(np.full_like(scatter[:,0], j+x_off), 0.03), scatter[:, nlab*i+j],
                                                s = 100, alpha = .5, c = col)

        ci = data.iloc[nlab*i:nlab*(i+1), 1:3].values.T
        errobars = data.iloc[nlab*i:nlab*(i+1), 3:].values.T
            
        plt.errorbar(x_loc, bar_height,
                    yerr = errobars, c = 'none', ecolor='k',
                    elinewidth=4, capsize = 10, capthick=4)

        if HIGHLIGHT:
            #is_sig = np.logical_and(bar_height > 0, np.abs(errobars[0]) < bar_height)
            #is_sig += np.logical_and(bar_height < 0, np.abs(errobars[1]) < -bar_height)
            is_sig = np.logical_and(bar_height > 0, ci[0] > 0)
            is_sig += np.logical_and(bar_height < 0, ci[1] < 0)
            is_sig = np.where(is_sig > 0)[0]
            #is_sig = np.arange(nlab)
            for j in [-1, 1]:
                plt.bar(x_loc[is_sig], j, width = width, color = 'yellow', alpha = .3, zorder = 0)
        
    plt.axhline(y = 0,linewidth = 4, color = 'k')

    plt.xticks(np.arange(nlab), labels)
    plt.tick_params(labelsize = 40)
    ax.tick_params(axis='x', pad=20)

    plt.setp(ax.get_xticklabels(), rotation=rot, ha="center", rotation_mode="anchor")

    #plt.xlabel('Behavioral Variables', fontsize = 60)
    plt.ylabel(ylabels, fontsize = 60)

    plt.legend(fontsize = 45, loc = 'best')

    plt.tight_layout()

    #if SAVE:
    #    plt.savefig('/Users/AlexCionca/Google Drive/My Drive/_Work/HUG/Data/Results/Severity_PLSC/'+
    #                '{}_Loading_V01.png'.format(study))
    return

def plot_behaviour(study, path_to_study, LC_id = 0, results = "loadings", invert = False, imaging = False,
                    labels = ['Obj. fatigue', 'Subj. fatigue', 'Cognitive', 'Physical', 'Social', 'Psychiatric', 'Age', 'Gender'],
                    groups = ['NEU-COV', 'ALL-COV', 'PSY-COV'],
                    colors = ['tab:red', 'grey', 'tab:blue'], figsize = (35, 20), **kwargs):

    data_type = 'Design'
    data_type_dim = 'behav'
    if imaging:
        data_type = 'Imaging'
        data_type_dim = 'img'

    norm_param = 2
    if op.isfile(op.join(path_to_study, study, 'myPLS_behavior_norm1-1_res.mat')):
        norm_param = 1

    with h5py.File(op.join(path_to_study, study, f'myPLS_behavior_norm{norm_param}-{norm_param}_res.mat'), 'r') as f:
        try:
            loadings_boot = f['res']['boot_results'][f'LC_{data_type_dim}_{results}_boot'][:,LC_id,:]
        except KeyError:
            print(f'Bootstraped {results} could not be found. Were they saved?')
            loadings_boot = None
        #behav_names = f['res']['design_names'][:]
        #groups = f['res']['group_names'][:]

    #print(behav_names, groups)

    res_w_cap = 'Loadings'
    if results == 'saliences':
        res_w_cap = 'Saliences'


    filename = f'myPLS_behavior_norm{norm_param}-{norm_param}_LC{LC_id+1}_{data_type}_{res_w_cap}.csv'
    cog_loadings = pd.read_csv(op.join(path_to_study, study, filename))

    #cog_loadings = pd.read_csv('/Users/AlexCionca/Documents/MATLAB/myPLS/{}/'.format(which_analysis)+
    #                           'myPLS_behavior_norm1-1_LC1_Imaging_Loadings.csv')

    #cog_loadings['Loadings'] = f['res']['LC_behav_loadings'][:][LC_id]
    #cog_loadings['Lower_CI'] = f['res']['boot_results']['LC_behav_loadings_lB'][:][LC_id]
    #cog_loadings['Upper_CI'] = f['res']['boot_results']['LC_behav_loadings_uB'][:][LC_id]

    cog_loadings['err_l'] = cog_loadings.Loadings - cog_loadings.Lower_CI
    cog_loadings['err_u'] = cog_loadings.Upper_CI - cog_loadings.Loadings

    # FATIGUE
    #labels = ['Obj. fatigue', 'Subj. fatigue', 'Cognitive', 'Physical', 'Social', 'Psychiatric', 'Age', 'Gender']
    #groups = ['NEU-COV', 'ALL-COV', 'PSY-COV']
    #colors = ['tab:red', 'grey', 'tab:blue']

    #labels = ['Obj. fatigue', 'Subj. fatigue', 'Cognitive', 'Physical', 'Social', 'Psychiatric', 'Age', 'Gender']
    #groups = ['NEU-COV', 'ALL-COV', 'PSY-COV']
    #colors = ['tab:red', 'grey', 'tab:blue']

    if invert:
        cog_loadings = -cog_loadings
        loadings_boot = -loadings_boot

    plot_bar(cog_loadings, labels, scatter = loadings_boot, groupings=groups, colors=colors,
                fig = None, ax = None, figsize = figsize, **kwargs)

    return


