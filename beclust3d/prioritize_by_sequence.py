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
from scipy.stats import norm

from _prioritize_by_sequence_plots_ import *

def prioritize_by_sequence(
    df_struc, df_consrv, df_nomutation, 
    workdir, 
    input_gene, input_screen, structureid, 
    function_type='mean', 
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
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        structureid: str, required
            the name of the AF and uniprot input

    Returns
        df_struc_consvr: pandas dataframe
            DataFrame
    """
    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir)
    edits_filedir = edits_filedir / input_gene
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')

    df_struc_consvr = df_struc
    if df_consrv is None: 
        df_struc_consvr['human_res_pos'] = df_struc_consvr['unipos']
        df_struc_consvr['conservation']  = 'None'
    else: 
        df_struc_consvr['human_res_pos'] = df_consrv['human_res_pos']
        df_struc_consvr['mouse_res_pos'] = df_consrv['mouse_res_pos']
        df_struc_consvr['mouse_res']     = df_consrv['mouse_res']
        df_struc_consvr['conservation']  = df_consrv['conservation']
    del df_struc, df_consrv

    # ALLOW USER TO PICK AGGREGATION FUNCTION #
    match function_type: 
        case 'mean': function = statistics.mean
        case 'median': function = statistics.median
        case 'min': function = min
        case 'max': function = max
        case _: warnings.warn('Warning: Invalid Function Type (mean, median, min, max)')

    struc_consrv_filename =  f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    df_struc_consvr.to_csv(edits_filedir / struc_consrv_filename, sep = "\t", index=False)

    # FOR EACH EDIT TYPE, AGGREGATE LFC AND EDITS WITH CONSERVATION #
    for edit_type in ['Missense', 'Silent', 'Nonsense']: 

        in_filename = f"screendata/{input_gene}_{screen_name}_{edit_type}_edits_list.tsv"
        df_edit = pd.read_csv(edits_filedir / in_filename, sep='\t')
        
        arr_unique_LFC = []
        arr_all_edits = []

        # FOR EACH RESIDUE #
        for i in range(len(df_struc_consvr)): 
            human_res_pos = df_struc_consvr.at[i, 'human_res_pos'] ### rewrite, and also rename human_res_pos
            df_pos_edits = df_edit.loc[df_edit['edit_pos'] == int(human_res_pos), ].reset_index() ### rewrite

            if len(df_pos_edits) > 1: 
                pos_LFCscore_list = df_pos_edits['LFC'].tolist()
                unique_LFC_res = round(function(pos_LFCscore_list), 3)

                pos_edits_list = df_pos_edits['this_edit'].tolist()
                all_edits_res = ';'.join(list(set(pos_edits_list)))
            elif len(df_pos_edits) == 1:   
                unique_LFC_res = round(df_pos_edits.at[0, 'LFC'], 3)
                all_edits_res = df_pos_edits.at[0, 'this_edit']
            else:
                unique_LFC_res, all_edits_res = '-', '-'

            arr_unique_LFC.append(unique_LFC_res)
            arr_all_edits.append(all_edits_res)

        df_struc_consvr[f'mean_{edit_type}_LFC'] = arr_unique_LFC
        df_struc_consvr[f'all_{edit_type}_edits'] = arr_all_edits

        # CALCULATE Z SCORE #

        list_z_LFC, list_p_LFC, list_plab_LFC = [], [], []
        # negative
        df_nomutation_neg = df_nomutation.loc[df_nomutation['LFC'] < 0.0, ]
        mu_neg = df_nomutation_neg['LFC'].mean()
        sigma_neg = df_nomutation_neg['LFC'].std()
        # positive
        df_nomutation_pos = df_nomutation.loc[df_nomutation['LFC'] > 0.0, ]
        mu_pos = df_nomutation_pos['LFC'].mean()
        sigma_pos = df_nomutation_pos['LFC'].std()

        for i in range(len(df_struc_consvr)):
            LFC_raw = df_struc_consvr.at[i, f'mean_{edit_type}_LFC']

            if LFC_raw == '-':
                LFC, z_LFC, p_LFC, plab_LFC = 0.0, '-', 1.0, 'p=1.0'
            else:
                LFC = float(df_struc_consvr.at[i, f'mean_{edit_type}_LFC'])

                if (LFC < 0.0):
                    z_LFC = statistics.NormalDist(mu=mu_neg, sigma=sigma_neg).zscore(LFC)
                    p_LFC = norm.sf(abs(z_LFC))
                    if            z_LFC < -3.29: plab_LFC = '-p=0.001'
                    elif -3.29 <= z_LFC < -2.58: plab_LFC = '-p=0.01'
                    elif -2.58 <= z_LFC < -1.96: plab_LFC = '-p=0.05'
                    elif -1.96 <= z_LFC < -1.65: plab_LFC = '-p=0.1'
                    elif -1.65 <= z_LFC < -1.0:  plab_LFC = '-p=0.3'
                    else:                        plab_LFC = '-p=1.0'

                elif (LFC > 0.0):
                    z_LFC = statistics.NormalDist(mu=mu_pos, sigma=sigma_pos).zscore(LFC)
                    p_LFC = norm.sf(abs(z_LFC))
                    if   3.29 < z_LFC:         plab_LFC = '+p=0.001'
                    elif 2.58 < z_LFC <= 3.29: plab_LFC = '+p=0.01'
                    elif 1.96 < z_LFC <= 2.58: plab_LFC = '+p=0.05'
                    elif 1.65 < z_LFC <= 1.95: plab_LFC = '+p=0.1'
                    elif  1.0 < z_LFC <= 1.65: plab_LFC = '+p=0.3'
                    else:                      plab_LFC = '+p=1.0'

                else: plab_LFC = 'p=1.0'

            list_z_LFC.append(z_LFC)
            list_p_LFC.append(p_LFC)
            list_plab_LFC.append(plab_LFC)

        df_struc_consvr[f'mean_{edit_type}_LFC_Z'] = list_z_LFC
        df_struc_consvr[f'mean_{edit_type}_LFC_p'] = list_p_LFC
        df_struc_consvr[f'mean_{edit_type}_LFC_plab'] = list_plab_LFC
        df_struc_consvr.round(4)

        # PLOT SCATTERPLOT AND COUNTS PLOT #
        if edit_type == 'Missense': 
            counts_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type, )
            scatterplot_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type, function_type, )
            scatterplot_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type, function_type, input='_Z')
            dual_scatterplot_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name)
            dual_histogram_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name)

    strcons_edits_filename = f"screendata/{input_gene}_{screen_name}_struc_consrv_proteinedits.tsv"
    df_struc_consvr.to_csv(edits_filedir / strcons_edits_filename, sep = '\t', index=False)

    return df_struc_consvr
