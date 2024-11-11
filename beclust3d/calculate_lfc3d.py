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
        workdir, input_gene, input_screens, screen_names=[], 
        nRandom=1000, function_type='mean', mut='Missense', 
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
        df_struct_3d: pandas dataframe
            a dataframe listing calculated LFC3D scores and their randomizations
    """

    edits_filedir = Path(workdir + '/' +  input_gene)
    if not screen_names: # for screen_names being an empty list
        screen_names = [input_screen.split('.')[0] for input_screen in input_screens]
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'LFC3D'):
        os.mkdir(edits_filedir / 'LFC3D')

    df_struct_3d = pd.DataFrame()
    df_struct_3d['unipos'] = df_str_cons['unipos']
    df_struct_3d['unires'] = df_str_cons['unires']

    # FOR EVERY SCREEN #
    for input_screen, screen_name in zip(input_screens, screen_names):

        str_cons_filename = edits_filedir / f"screendata/{input_gene}_{screen_name}_proteinedits.tsv"
        if not os.path.exists(str_cons_filename): 
            warnings.warn(f"{str_cons_filename} does not exist")
            if len(input_screens) == 1: 
                return None
            continue
        df_struc_edits = pd.read_csv(str_cons_filename, sep = "\t")

        taa_wise_norm_LFC = []
        for aa in range(len(df_struc_edits)): # for every residue
            taa_naa_wBE_LFC, sum_taa_naa_LFC = helper(df_struc_edits, aa, lookup=f'{function_type}_{mut}_LFC')
            if taa_naa_wBE_LFC == 0:
                taa_wise_norm_LFC.append('-')
            else: 
                taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))
        
        df_struct_3d[f"{screen_name}_LFC"] = df_struc_edits[f'{function_type}_{mut}_LFC']
        df_struct_3d[f"{screen_name}_LFC_Z"] = df_struc_edits[f'{function_type}_{mut}_LFC_Z']
        df_struct_3d[f"{screen_name}_LFC3D"] = taa_wise_norm_LFC
    
        str_cons_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name}_{mut}_proteinedits_rand.tsv"
        screen_name = input_screen.split('.')[0]
        if not os.path.exists(str_cons_filename): 
            warnings.warn(f"{str_cons_filename} does not exist")
            if len(input_screens) == 1: 
                return None
            continue
        
        df_struc_edits = pd.read_csv(str_cons_filename, sep = "\t")
        dict_temp = {}
        for r in range(0, nRandom):
            taa_wise_norm_LFC = []

            for aa in range(0, len(df_struc_edits)):
                taa_naa_wBE_LFC, sum_taa_naa_LFC = helper(df_struc_edits, aa, lookup=f'{function_type}_missense_LFC')
                if taa_naa_wBE_LFC == 0:
                    taa_wise_norm_LFC.append('-')
                else: 
                    taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))

            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = df_struc_edits[f'{function_type}_missense_LFCr{str(r+1)}']
            dict_temp[f"{screen_name}_LFC3Dr{str(r+1)}"] = taa_wise_norm_LFC

        df_struct_3d = pd.concat((df_struct_3d, pd.DataFrame(dict_temp)), axis=1)

    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr.tsv"
    df_struct_3d.to_csv(out_filename, sep = '\t', index=False)

    return df_struct_3d

def helper(
    df_struc_edits, aa, lookup
): 
    taa_naa_wBE_LFC, sum_taa_naa_LFC = 0, 0.0

    naa_pos_str = df_struc_edits.at[aa, 'Naa_pos']
    # there are no residues nearby for low radius
    if not isinstance(naa_pos_str, float): 
        naa_pos_list = naa_pos_str.split(';') # neighboring residue positions
        taa_LFC = df_struc_edits.at[aa, lookup] # target LFC

        if taa_LFC != '-':
            taa_naa_wBE_LFC, sum_taa_naa_LFC = 1, float(taa_LFC)

        for naa_pos in naa_pos_list: 
            naa_LFC = df_struc_edits.at[int(naa_pos)-1, lookup]
            if naa_LFC != '-': 
                sum_taa_naa_LFC += float(naa_LFC)
                taa_naa_wBE_LFC += 1
    else: 
        taa_LFC = df_struc_edits.at[aa, lookup] # target LFC
        if taa_LFC != '-':
            taa_naa_wBE_LFC, sum_taa_naa_LFC = 1, float(taa_LFC)

    return taa_naa_wBE_LFC, sum_taa_naa_LFC
