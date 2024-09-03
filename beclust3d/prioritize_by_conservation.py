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
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def prioritize_by_conservation(
        df_struc, df_consrv, 
        workdir, 
        input_gene, input_screen, structureid, 
        function_type='mean', 
): 
    """
    Description
        Takes in results across multiple edit types for a screen, and
        aggregates the edits for each residue with conservation information. 

    Params
        df_struc: pandas dataframe, required
            DataFrame output from conservation()
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
                unique_LFC_res = str(round(function(pos_LFCscore_list), 3))

                pos_edits_list = df_pos_edits['this_edit'].tolist()
                all_edits_res = ';'.join(list(set(pos_edits_list)))
            elif len(df_pos_edits) == 1:   
                unique_LFC_res = str(round(df_pos_edits.at[0, 'LFC'], 3))
                all_edits_res = df_pos_edits.at[0, 'this_edit']
            else:
                unique_LFC_res, all_edits_res = '-', '-'

            arr_unique_LFC.append(unique_LFC_res)
            arr_all_edits.append(all_edits_res)

        df_struc_consvr[f'mean_{edit_type}_LFC'] = arr_unique_LFC
        df_struc_consvr[f'all_{edit_type}_edits'] = arr_all_edits

        # PLOT SCATTERPLOT AND COUNTS PLOT #
        if edit_type == 'Missense': 
            counts_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type, )
            scatterplot_by_residue(df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type, function_type, )

    strcons_edits_filename = f"screendata/{input_gene}_{screen_name}_struc_consrv_proteinedits.tsv"
    df_struc_consvr.to_csv(edits_filedir / strcons_edits_filename, sep = '\t', index=False)

    return df_struc_consvr


def counts_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, 
    edit_type, 
): 
    # PREP DATA #
    counts = df_struc_consvr[f'all_{edit_type}_edits'].str.count(';').fillna(0).astype(int)+1
    counts[df_struc_consvr[f'all_{edit_type}_edits'] == '-'] = 0

    # PLOT #
    plt.figure(figsize=(10, 4))
    ax = sns.barplot(x=df_struc_consvr['unipos'], y=counts, 
                     color='steelblue', edgecolor='steelblue')
    ax.set_ylabel(f"Count of {edit_type} Mutations")
    ax.set_title(f"Count of {edit_type} Mutations Per Residue {screen_name}")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    counts_filename = f"plots/{input_gene}_{screen_name}_num_{edit_type}_per_residue.pdf"
    plt.savefig(edits_filedir / counts_filename, dpi=300)

def scatterplot_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, 
    edit_type, function_type, 
): 
    # PREP DATA #
    x_list = df_struc_consvr['unipos'].tolist()
    y_list = df_struc_consvr[f'{function_type}_{edit_type}_LFC'].tolist()
    x_vals = [x for x, y in zip(x_list, y_list) if y!='-']
    y_vals = [float(y) for y in y_list if y!='-']

    # PLOT #
    plt.figure(figsize=(10, 4))
    ax = sns.scatterplot(x=x_vals, y=y_vals, color='steelblue', edgecolor='steelblue')
    ax.axhline(-1.0, c="red", linestyle="--")
    ax.axhline(1.0, c="blue", linestyle="--")
    ax.axhline(0.0, c="gray", linestyle="--")
    ax.set_ylabel(f"{edit_type} LFC Score")
    ax.set_title(f'{edit_type} LFC Score By Residue {screen_name}')
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    scatter_filename = f"plots/{input_gene}_{screen_name}_{edit_type}_lfc_score_by_residue.pdf"
    plt.savefig(edits_filedir / scatter_filename, dpi=300)
