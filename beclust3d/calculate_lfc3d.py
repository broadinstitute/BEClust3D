"""
File: calculate_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3

"""

import pandas as pd
from pathlib import Path
import os
import warnings

def calculate_lfc3d(
        df_str_cons, 
        workdir, input_gene, input_screens, 
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
        input_screens: str, required
            the name of the input screen
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_str_cons_3daggr: pandas dataframe
            a dataframe listing calculated LFC3D scores and their randomizations
    """

    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'LFC3D'):
        os.mkdir(edits_filedir / 'LFC3D')

    df_str_cons_3daggr = pd.DataFrame()
    df_str_cons_3daggr['unipos'] = df_str_cons['unipos']
    df_str_cons_3daggr['unires'] = df_str_cons['unires']

    for input_screen in input_screens: # for every screen

        screen_name = input_screen.split('.')[0]
        str_cons_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name}_struc_consrv_missenseedits_randomized.tsv"
        if not os.path.exists(str_cons_filename): 
            warnings.warn(f"{str_cons_filename} does not exist")
            if len(input_screens) == 1: 
                return None
            continue
        df_str_cons_edits = pd.read_csv(str_cons_filename, sep = "\t")

        taa_wise_norm_LFC = []

        for aa in range(0, len(df_str_cons_edits)): # for every residue
            taa_naa_wBE_LFC, sum_taa_naa_LFC = helper(df_str_cons_edits, aa)
            if taa_naa_wBE_LFC == 0:
                taa_wise_norm_LFC.append('-')
            else: 
                taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))
        
        df_str_cons_3daggr[f"{screen_name}_LFC"] = df_str_cons_edits['mean_Missense_LFC']
        df_str_cons_3daggr[f"{screen_name}_LFC_Z"] = df_str_cons_edits['mean_Missense_LFC_Z']
        df_str_cons_3daggr[f"{screen_name}_LFC3D"] = taa_wise_norm_LFC
    
        dict_temp = {}
        for r in range(0, nRandom):
            taa_wise_norm_LFC = []

            for aa in range(0, len(df_str_cons_edits)):
                taa_naa_wBE_LFC, sum_taa_naa_LFC = helper(df_str_cons_edits, aa)
                if taa_naa_wBE_LFC == 0:
                    taa_wise_norm_LFC.append('-')
                else: 
                    taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))

            ### i think this randomization lookup just pull from previous shuffled data, so you don't need to call the dataframe but rather just copy over columns
            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = df_str_cons_edits[f'mean_missense_LFCr{str(r+1)}'] ###
            dict_temp[f"{screen_name}_LFC3Dr{str(r+1)}"] = taa_wise_norm_LFC

        df_str_cons_3daggr = pd.concat((df_str_cons_3daggr, pd.DataFrame(dict_temp)), axis=1)

    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_per_Random_LFC3Dr.tsv"
    df_str_cons_3daggr.to_csv(out_filename, sep = '\t', index=False)

    return df_str_cons_3daggr

def helper(
    df_str_cons_edits, aa
): 
    taa_naa_wBE_LFC = 0
    sum_taa_naa_LFC = 0.0

    naa_pos_str = df_str_cons_edits.at[aa, 'Naa_pos']
    # there are no residues nearby for low radius
    if not isinstance(naa_pos_str, float): 
        naa_pos_list = naa_pos_str.split(';') # neighboring residue positions
        taa_LFC = df_str_cons_edits.at[aa, 'mean_missense_LFC'] # target LFC

        if taa_LFC != '-':
            taa_naa_wBE_LFC = 1
            sum_taa_naa_LFC = float(taa_LFC)

        for naa_pos in naa_pos_list: 
            naa_LFC = df_str_cons_edits.at[int(naa_pos)-1, 'mean_missense_LFC']
            if naa_LFC != '-': 
                sum_taa_naa_LFC += float(naa_LFC)
                taa_naa_wBE_LFC += 1
    else: 
        taa_LFC = df_str_cons_edits.at[aa, 'mean_missense_LFC'] # target LFC
        if taa_LFC != '-':
            taa_naa_wBE_LFC = 1
            sum_taa_naa_LFC = float(taa_LFC)

    return taa_naa_wBE_LFC, sum_taa_naa_LFC
