"""
File: calculate_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3

"""

import pandas as pd
from pathlib import Path
import os

def calculate_lfc3d(
        df_str_cons, 
        workdir, 
        input_gene, input_screen, 
        nRandom=1000, 
): 
    """
    Description
        Calculates LFC 3D scores from structural conservation data

    Params
        df_str_cons: pandas dataframe
            from previous step randomize_by_conservation()
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_str_cons_3daggr: pandas dataframe
            a dataframe listing calculated LFC3D scores and their randomizations
    """

    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'LFC3D'):
        os.mkdir(edits_filedir / 'LFC3D')
    
    str_cons_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name}_struc_consrv_missenseedits_randomized.tsv"
    df_str_cons_edits = pd.read_csv(str_cons_filename, sep = "\t")
    
    df_str_cons_3daggr = pd.DataFrame()
    df_str_cons_3daggr['unipos'] = df_str_cons['unipos']
    df_str_cons_3daggr['unires'] = df_str_cons['unires']
    taa_wise_norm_LFC = []

    for aa in range(0, len(df_str_cons_edits)):
        taa_naa_wBE_LFC, sum_taa_naa_LFC = helper(df_str_cons_edits, aa)
        if taa_naa_wBE_LFC == 0:
            taa_wise_norm_LFC.append('-')
        else:
            taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))
    
    df_str_cons_3daggr[f"{screen_name}_LFC"] = df_str_cons_edits['mean_missense_LFC']
    df_str_cons_3daggr[f"{screen_name}_LFC3D"] = taa_wise_norm_LFC
    
    dict_temp = {}
    for r in range(0, nRandom):
        col_head = 'mean_missense_LFCr' + str(r+1)
        taa_wise_norm_LFC = []

        for aa in range(0, len(df_str_cons_edits)):
            taa_naa_wBE_LFC, sum_taa_naa_LFC = helper(df_str_cons_edits, aa)
            if taa_naa_wBE_LFC == 0:
                taa_wise_norm_LFC.append('-')
            else:
                taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))

        dict_temp[f"{screen_name}_LFC{str(r+1)}"] = df_str_cons_edits[col_head]
        dict_temp[f"{screen_name}_LFC3D{str(r+1)}"] = taa_wise_norm_LFC

    df_str_cons_3daggr = pd.concat((df_str_cons_3daggr, pd.DataFrame(dict_temp)), axis=1)
    out_filename = edits_filedir / f"LFC3D/{input_gene}_{screen_name}_LFC_LFC3D_per_Random_LFC3Dr.tsv"
    df_str_cons_3daggr.to_csv(out_filename, sep = '\t', index=False)

    return df_str_cons_3daggr

def helper(
    df_str_cons_edits, aa
): 
    naa_list = df_str_cons_edits.at[aa, 'Naa'].split(';') # neighboring amino acids
    naa_pos_list = df_str_cons_edits.at[aa, 'Naa_pos'].split(';') # neighboring residue positions
    taa_LFC = df_str_cons_edits.at[aa, 'mean_missense_LFC'] # target LFC

    taa_naa_wBE_LFC = 0 # residues that are conserved and with a BE edit value
    sum_taa_naa_LFC = 0.0
    if taa_LFC != '-':
        taa_naa_wBE_LFC = 1
        sum_taa_naa_LFC = float(taa_LFC)

    for j in range(0, len(naa_list)):
        naa_pos = int(naa_pos_list[j])
        naa_LFC = df_str_cons_edits.at[naa_pos-1, 'mean_missense_LFC']
        if naa_LFC != '-':
            sum_taa_naa_LFC = sum_taa_naa_LFC + float(naa_LFC)
            taa_naa_wBE_LFC = taa_naa_wBE_LFC + 1

    return taa_naa_wBE_LFC, sum_taa_naa_LFC
