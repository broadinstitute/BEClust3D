"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2

"""

import pandas as pd
from pathlib import Path
import os
import statistics
import warnings
import numpy as np
from scipy.stats import norm
warnings.filterwarnings("ignore", category=RuntimeWarning)

from _prioritize_by_sequence_plots_ import *


def get_plabel(z_LFC, direction):
    if direction == 'negative': 
        thresholds = [(-3.29,'-p=0.001'), (-2.58,'-p=0.01'), (-1.96,'-p=0.05'), (-1.65,'-p=0.1'), (-1.0,'-p=0.3')] 
    else: 
        thresholds = [(3.29,'+p=0.001'), (2.58,'+p=0.01'), (1.96,'+p=0.05'), (1.65,'+p=0.1'), (1.0,'+p=0.3')]
    for threshold, label in thresholds:
        if (direction == 'negative' and z_LFC < threshold) or (direction == 'positive' and z_LFC > threshold):
            return label
    return '-p=1.0' if direction=='negative' else '+p=1.0'

def prioritize_by_sequence(
    df_struc, df_consrv, df_nomutation, 
    workdir, 
    input_gene, screen_name, file_dict, 
    function=statistics.mean, function_name='mean', target_res_pos='human_res_pos', 
): 
    """
    Description
        Takes in results across multiple edit types for a screen, and
        aggregates the edits for each residue with sequence and conservation information. 

    Params
        df_struc: pandas dataframe, required
            DataFrame output from af_structural_features()
        df_consrv: pandas dataframe, required
            DataFrame output from conservation()
        df_nomutation: pandas dataframe, required
            DataFrame output from preprocess_be_results() to calculate baseline on
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        structureid: str, required
            the name of the AF and uniprot input
        screen_name: str, optional
            the name of the input screen
        function: function, optional
            indicates the type of supported aggregation function for a residue
            options: statistics.mean, statistics.median, min, max

    Returns
        df_protein: pandas dataframe
            DataFrame
    """
    
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')

    df_protein = df_struc
    if df_consrv is None: 
        df_protein['human_res_pos'] = df_protein['unipos']
        df_protein['conservation']  = 'None'
    else: 
        df_protein['human_res_pos'] = df_consrv['human_res_pos']
        df_protein['mouse_res_pos'] = df_consrv['mouse_res_pos']
        df_protein['mouse_res']     = df_consrv['mouse_res']
        df_protein['conservation']  = df_consrv['conservation']
    del df_struc, df_consrv

    # struc_consrv_filename =  f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    # df_protein.to_csv(edits_filedir / struc_consrv_filename, sep = "\t", index=False)

    # FOR EACH EDIT TYPE, AGGREGATE LFC AND EDITS WITH CONSERVATION #
    for in_filename, mut in file_dict.items(): 
        df_edit = pd.read_csv(edits_filedir / in_filename, sep='\t')
        
        arr_unique_LFC = []
        arr_unique_LFC_stdev = []
        arr_all_edits = []

        # FOR EACH RESIDUE #
        for i in range(len(df_protein)): 
            human_res_pos = df_protein.at[i, target_res_pos]
            df_pos_edits = df_edit.loc[df_edit['edit_pos'] == int(human_res_pos), ].reset_index() ###

            ### should try to rework logic here 250205
            if (df_consrv is None) or (df_protein.at[i, 'mouse_res'] != '-'): 
                if len(df_pos_edits) > 1: 
                    score_list = df_pos_edits['LFC'].tolist()
                    unique_LFC_res = round(function(score_list), 3)
                    stdev_res = np.std(score_list)

                    pos_edits_list = df_pos_edits['this_edit'].tolist()
                    all_edits_res = ';'.join(list(set(pos_edits_list)))
                elif len(df_pos_edits) == 1: 
                    unique_LFC_res = round(df_pos_edits.at[0, 'LFC'], 3)
                    stdev_res = 0
                    all_edits_res = df_pos_edits.at[0, 'this_edit']
                else:
                    unique_LFC_res, all_edits_res, stdev_res = '-', '-', '-'
            else:
                unique_LFC_res, all_edits_res, stdev_res = '-', '-', '-'

            arr_unique_LFC.append(unique_LFC_res)
            arr_unique_LFC_stdev.append(stdev_res)
            arr_all_edits.append(all_edits_res)

        df_protein[f'{function_name}_{mut}_LFC'] = arr_unique_LFC
        df_protein[f'{function_name}_{mut}_LFC_stdev'] = arr_unique_LFC_stdev
        df_protein[f'all_{mut}_edits'] = arr_all_edits

        # CALCULATE Z SCORE #
        # FOR NEG AND POS SEPARATELY, CALC Z SCORE BASED ON THE MEAN STD PER SCREEN PER GENE PER DIRECTION #
        neg_mask = df_nomutation['LFC'] < 0.0 # NEG #
        pos_mask = df_nomutation['LFC'] > 0.0 # POS #
        mu_neg, sigma_neg = df_nomutation.loc[neg_mask, 'LFC'].mean(), df_nomutation.loc[neg_mask, 'LFC'].std()
        mu_pos, sigma_pos = df_nomutation.loc[pos_mask, 'LFC'].mean(), df_nomutation.loc[pos_mask, 'LFC'].std()

        list_z_LFC, list_p_LFC, list_plab_LFC = [], [], []

        for i in range(len(df_protein)):
            LFC_raw = df_protein.at[i, f'{function_name}_{mut}_LFC']

            if LFC_raw == '-': 
                LFC, z_LFC, p_LFC, plab_LFC = 0.0, '-', 1.0, 'p=1.0'
            else: 
                LFC = float(df_protein.at[i, f'{function_name}_{mut}_LFC'])

                if (LFC < 0.0):
                    z_LFC = statistics.NormalDist(mu=mu_neg, sigma=sigma_neg).zscore(LFC)
                    p_LFC = norm.sf(abs(z_LFC))
                    plab_LFC = get_plabel(z_LFC, direction='negative')
                elif (LFC > 0.0):
                    z_LFC = statistics.NormalDist(mu=mu_pos, sigma=sigma_pos).zscore(LFC)
                    p_LFC = norm.sf(abs(z_LFC))
                    plab_LFC = get_plabel(z_LFC, direction='positive')
                else: plab_LFC = 'p=1.0' ### redundant? 250121

            list_z_LFC.append(z_LFC)
            list_p_LFC.append(p_LFC)
            list_plab_LFC.append(plab_LFC)

        df_protein[f'{function_name}_{mut}_LFC_Z'] = list_z_LFC
        df_protein[f'{function_name}_{mut}_LFC_p'] = list_p_LFC
        df_protein[f'{function_name}_{mut}_LFC_plab'] = list_plab_LFC
        df_protein.round(4)

    strcons_edits_filename = f"screendata/{input_gene}_{screen_name}_proteinedits.tsv"
    df_protein.to_csv(edits_filedir / strcons_edits_filename, sep = '\t', index=False)

    return df_protein

def plots_by_sequence(
    df_protein, 
    workdir, 
    input_gene, screen_name, function_name='mean', 
): 
    edits_filedir = Path(workdir)
    # PLOT SCATTERPLOT AND COUNTS PLOT #
    counts_by_residue(df_protein, edits_filedir, input_gene, screen_name, 'Missense', )
    stdev_by_residue(df_protein, edits_filedir, input_gene, screen_name, function_name, 'Missense')
    stdev_by_residue(df_protein, edits_filedir, input_gene, screen_name, function_name, 'Missense', yaxis=False)
    scatterplot_by_residue(df_protein, edits_filedir, input_gene, screen_name, 'Missense', function_name, )
    scatterplot_by_residue(df_protein, edits_filedir, input_gene, screen_name, 'Missense', function_name, input='_Z')
    dual_scatterplot_by_residue(df_protein, edits_filedir, input_gene, screen_name, function_name)
    dual_histogram_by_residue(df_protein, edits_filedir, input_gene, screen_name, function_name)
