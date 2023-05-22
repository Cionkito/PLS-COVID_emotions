#!/usr/bin/env python
# coding: utf-8
""" Python module with utility functions for handling and analysing fMRI data

Run as:

    python Visualize/plot_PLSC.py
"""

import numpy as np
import pandas as pd

import os.path as op

def get_roi_name(labelname):
    """ Get the ROI name following the "Hemisphere_Network_*Region*_Suffix" pattern.

    Parameters
    ----------
    labelname : str
        Name of the label

    Returns
    -------
    str
        Name of the ROI extracted from the label
    """
    return str.join('_', labelname.split('_')[2:])

def fix_Vermis(roi_name, network):
    """ Fix the ROI name and add network label

    Parameters
    ----------
    roi_name : str
        Name of the ROI to fix
    network : str
        Name of the network to append

    Returns
    -------
    str
        Fixed name of the ROI with added network
    """

    if "Vermis" in roi_name:
        fixed_name = str.join("_", ["", network, roi_name])
    else:
        splitted_name = roi_name.split("_")
        fixed_name = str.join("_", [splitted_name[0], network]+splitted_name[1:])

    return fixed_name

def get_atlas_labels(path_to_atlas, filters = "all", atlas_fname = "Yeo_WholeBrain_156_Labels.csv"):
    """ Get label of ROI from a custom atlas

    Parameters
    ----------
    path_to_atlas : str
        Path to the atlas directory
    filters : str, optional
        Filter to consider a subset of ROIs, by default "all"
    atlas_fname : str, optional
        Name of file with atlas labels

    Returns
    -------
    Pandas DataFrame
        Index, labels and anatomical locations of the atlas ROIs

    Raises
    ------
    FileNotFoundError
        Label file was not found.
    """

    if type(filters) != list:
        filters = [filters]

    n_reg = 0

    if "all" in filters:
        n_reg = 156
    else:
        n_reg = 156
        #n_reg = 100*("cortic" in fltr) + 34*("cereb" in fltr) + 22*("BG" in fltr)

    if 'csv' in atlas_fname:
        atlas_labels = pd.read_csv(op.join(op.realpath(path_to_atlas), atlas_fname), index_col = 0)
    if 'txt' in atlas_fname:
        atlas_labels = pd.read_csv(op.join(op.realpath(path_to_atlas), atlas_fname), index_col = 0, sep=' ', names=['roi_name', 'id'])

    if 'Yeo_WholeBrain_156_Labels' in atlas_fname:
        atlas_labels.rename(columns={"Network":"roi_name"}, inplace=True)

        # Fix Vermis labels not containing a LH or RH for hemisphere localisation and add network labels
        atlas_labels.iloc[100:134, 0] = atlas_labels.iloc[100:134, 0].apply(fix_Vermis, network = "Cereb")
        atlas_labels.iloc[134:156, 0] = atlas_labels.iloc[134:156, 0].apply(fix_Vermis, network = "Subcortical")

        atlas_labels["Network"] = atlas_labels.roi_name.apply(lambda x: x.split("_")[1])

    if 'JHU' in atlas_fname:
        atlas_labels.rename(columns={"Network":"roi_name"}, inplace=True)
        #atlas_labels = atlas_labels.drop(columns = ['Unnamed: 0'])

    return atlas_labels

def sort_atlas_labels(labels):
    """ Sort atlas label by network

    Parameters
    ----------
    labels : Pandas DataFrame
        Atlas labels

    Returns
    -------
    sorted_labels : Pandas DataFrame
        Sorted label dataframe
    sorting_indices : Numpy array
        Indices for future sortings
    """

    unique_networks = pd.Series(np.arange(len(labels.Network.unique()))+1, index = labels.Network.unique())

    sorting_indices = np.argsort(labels.Network.apply(lambda x: unique_networks[x]).values)

    sorted_labels = labels.iloc[sorting_indices]
    net_idx = sorted_labels.loc[:, "Network"].apply(lambda x: unique_networks[x]).values
    #sorted_labels.loc[:, "Network_id"] = net_idx
    sorted_labels.insert(len(sorted_labels.columns), "Network_id", net_idx)

    return sorted_labels, sorting_indices

def fix_covidcog_names(dataframe, location = 'CovidCog'):

    pat_map = {'P088':'R3603', 'P092':'R3618', 'P124':'R3674'}

    dataframe_updated = dataframe.copy()

    for i, row in dataframe.iterrows():
        if row[location] in ['P088', 'P092', 'P124']:
            dataframe_updated.loc[i, 'CovidCog_mri'] = pat_map[row[location]]

    return dataframe_updated

def get_subject_DF(n_subjects, path_to_subject_infos):

    subj_DF = pd.read_excel(path_to_subject_infos, sheet_name = 'Visual ToPy', usecols = 'A:P', nrows = 80)

    subj_DF.loc[:, 'CovidCog_mri'] = subj_DF.loc[:, 'CovidCog']
    subj_DF = fix_covidcog_names(subj_DF)
    #pat_map = {'P088':'R3603', 'P092':'R3618', 'P124':'R3674'}
    #for i, row in subj_DF.iterrows():
    #    if row.CovidCog in ['P088', 'P092', 'P124']:
    #        subj_DF.loc[i, 'CovidCog_mri'] = pat_map[row.CovidCog]

    #to_keep = subj_DF.EMO == 1
    to_keep = subj_DF.MRI == 1
    subj_DF = subj_DF.loc[to_keep].drop(columns = ['Patient ID', 'SEV'])\
                        .sort_values(by = 'CovidCog_mri')

    if n_subjects == 48:
        res_DF = subj_DF.loc[subj_DF.CovidCog != 'P070']
        #res_DF = subj_DF.loc[subj_DF.CovidCog != 'P199']
        res_DF = res_DF[res_DF.CovidCog != 'P054']
        res_DF = res_DF[res_DF['Phase 2'] >= 0].reset_index(drop = True)
    elif n_subjects == 49:
        res_DF = subj_DF.loc[subj_DF.CovidCog != 'P070']
        #res_DF = subj_DF.loc[subj_DF.CovidCog != 'P199']
        #res_DF = res_DF[res_DF.CovidCog != 'P054']
        res_DF = res_DF[res_DF['Phase 2'] >= 0].reset_index(drop = True)
    elif n_subjects == 76:
        res_DF = subj_DF.loc[[pat_id not in ['P070'] for pat_id in subj_DF.CovidCog]]
        
        is_p2 = res_DF[res_DF['Phase 2'] == 1]
        res_DF.loc[res_DF['Phase 2'] == 1, 'Phase 2'] = 0
        is_p2.CovidCog = is_p2.CovidCog.apply(lambda x: x+"P2")
        
        res_DF = pd.concat([res_DF, is_p2]).sort_values(by = ['CovidCog_mri', 'Phase 2']).reset_index(drop = True)
        res_DF['Phase 2'] += 2
    else:
        res_DF = subj_DF.copy()

    return res_DF

def get_behavioural_data(path_to_behaviour, **kwargs):
    """_summary_

    Parameters
    ----------
    path_to_behaviour : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """

    behav_DF = pd.read_excel(path_to_behaviour, na_values=['ND', 'NA', 'NS'], **kwargs)

    if 'cov' in behav_DF.iloc[1, 0].lower():
        print('Generating CovidCog column ...')
        behav_DF['CovidCog']=behav_DF.iloc[:, 0].apply(lambda x: x.split('_')[1])
    elif 'id' in behav_DF.columns[0].lower():
        print('No CovidCog column found !')
        behav_DF['CovidCog']=behav_DF.iloc[:, 0]

    return behav_DF.sort_values(by = ['CovidCog']).reset_index(drop = True)

def get_olfaction_data(path_to_olfaction_data, olfaction_filename='20211019_Raw neuropsychological+Psychiatric+GERT_COVIDCOG Project_V1+V2_IJA.xlsx',
                       threshold_hyposmic=12, threshold_anosmic=8):
    
    olfaction_DF=get_behavioural_data(op.join(path_to_olfaction_data, olfaction_filename),
                                      sheet_name='Raw total 6 and 12 months').loc[:, ['CovidCog', 'Total Anosmie']]
    olfaction_DF['Anosmia']=2
    olfaction_DF.loc[olfaction_DF.loc[:, 'Total Anosmie'] > threshold_hyposmic, 'Anosmia']=3
    olfaction_DF.loc[olfaction_DF.loc[:, 'Total Anosmie'] < threshold_anosmic, 'Anosmia']=1
    
    olfaction_DF['is_hyposmic']=1 + (olfaction_DF.loc[:, 'Anosmia'] > 2).astype(int)

    return olfaction_DF

def merge_DF_w_MRI(subj_df, mri_subj_list, fix_p2 = False, mri_names = False):

    mri_subj_list_fixed = mri_subj_list.copy()
    if 'sub' in mri_subj_list[0]:
        print('Warning remove "sub" from participant label !')
        mri_subj_list_fixed = [sub.split('sub-')[-1] for sub in mri_subj_list]
    if type(mri_subj_list) != pd.DataFrame:
        mri_df = pd.DataFrame(mri_subj_list_fixed, columns = ['CovidCog_mri'])
    else:
        mri_df = mri_subj_list_fixed.copy()

    merged_df = mri_df.merge(subj_df, left_on = 'CovidCog_mri', right_on = 'CovidCog'+'_mri'*mri_names, how = 'left')
    
    if fix_p2:
        merged_df['is_p2'] = np.zeros(len(merged_df))
        sub_with_P2 = [sub for sub in mri_subj_list_fixed if 'P2' in sub]
        for sub in sub_with_P2:
            data_to_fill = merged_df.loc[merged_df.CovidCog == sub[:-2]].iloc[:, 2:].values.tolist()
            if len(data_to_fill) == 0:
                data_to_fill = [subj_df.loc[subj_df.CovidCog == sub[:-2]].iloc[:, 1:].values.tolist()[0] + [0]]
            
            data_to_fill = [sub] + data_to_fill[0][:-1] + [1]
            merged_df.loc[merged_df.loc[:, 'CovidCog_mri'+'_x'*(mri_names)] == sub, merged_df.columns[1:]] = data_to_fill
            #merged_df.loc[merged_df.loc[:, 'CovidCog_mri'+'_x'*(not mri_names)] == sub, merged_df.columns[1:]] = data_to_fill

    return merged_df

def add_phenotype(subj_df, path_to_phenotype, left_on = 'CovidCog'):

    phenotype_df = pd.read_csv(path_to_phenotype)

    merged_df = subj_df.merge(phenotype_df, left_on = left_on, right_on = 'CovidCog', how = 'left')

    merged_df = merged_df.drop(columns = merged_df.loc[:, ['Unnamed' in col for col in merged_df.columns]].columns)

    return merged_df

def fix_region_names(roi):
    if type(roi) == list:
        fixed_names = [fix_region_names(reg) for reg in roi]
    else:
        no_hemi = roi.split('H_')[-1]
        # fix name to the format: Network_{Region1 Region2 ...}
        network_lab = no_hemi.split('_')[0]
        region_lab = '\_'.join(no_hemi.split('_')[1:])

        fixed_names = network_lab + '_{' + region_lab + '}'

    return fixed_names

def save_node_edge_files(func_con_matrix, path_to_atlas, study = '', saveloc = None, **kwargs):

    atlas_labels = get_atlas_labels(path_to_atlas=path_to_atlas)

    binary_matrix = (func_con_matrix != 0)
    matrix_mask = binary_matrix.any(axis = 1)

    roi = atlas_labels.loc[matrix_mask].copy()
    # Try using .loc[row_indexer,col_indexer] = value instead
    fixed_roi = roi.roi_name.apply(fix_region_names).tolist()
    roi.loc[:, 'roi_name'] = fixed_roi
    #colors	degree	reg_lab
    roi.loc[:, 'degree'] = binary_matrix[matrix_mask][:, matrix_mask].sum(axis = 1)
    unique_networks = roi.Network.unique()
    network_to_id = pd.Series(np.arange(len(unique_networks)), index = unique_networks)
    roi.loc[:, 'colors'] = roi.Network.apply(lambda x: network_to_id[x]+1).values

    # R	A	S	colors	degree	reg_lab
    save_node = roi.loc[:, ['R', 'A', 'S', 'colors', 'degree', 'roi_name']]
    save_edge = func_con_matrix[matrix_mask][:, matrix_mask]

    print(save_edge.shape)

    if saveloc is not None:
        save_node.to_csv(op.join(saveloc,f'node-{study}.node'), header=False, index=False, sep='\t')
        np.savetxt(op.join(saveloc,f'edge-{study}.edge'), save_edge)#, fmt = '%1.4f')

    return save_node