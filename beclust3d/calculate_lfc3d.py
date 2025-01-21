"""
File: calculate_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3

"""

import pandas as pd
import numpy as np
from pathlib import Path
import os
import warnings
import numpy as np
import warnings
warnings.filterwarnings('ignore')

def calculate_lfc3d(
        df_str_cons, 
        workdir, input_gene, screen_names, str_cons_filenames, str_cons_rand_filenames, 
        nRandom=1000, function_type='mean', mut='Missense', function_3Daggr=np.mean, LFC_only=False, 
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
        screen_names: list of str, required
            the names of the input screens
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_struct_3d: pandas dataframe
            a dataframe listing calculated LFC3D scores and their randomizations
    """

    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'LFC'):
        os.mkdir(edits_filedir / 'LFC')
    if not os.path.exists(edits_filedir / 'LFC3D') and not LFC_only:
        os.mkdir(edits_filedir / 'LFC3D')

    df_struct_3d = pd.DataFrame()
    df_struct_3d['unipos'] = df_str_cons['unipos']
    df_struct_3d['unires'] = df_str_cons['unires']

    # FOR EVERY SCREEN #
    for screen_name, filename, rand_filename in zip(screen_names, str_cons_filenames, str_cons_rand_filenames):
        print(screen_name)
        
        # GET LFC and LFC3D VALUES #
        if not os.path.exists(edits_filedir / filename): 
            warnings.warn(f"{filename} does not exist")
        df_struc_edits = pd.read_csv(edits_filedir / filename, sep = "\t")

        if not LFC_only: 
            taa_wise_norm_LFC = []
            for aa in range(len(df_struc_edits)): # FOR EVERY RESIDUE #
                taa_naa_LFC_vals = helper(df_struc_edits, aa, lookup=f'{function_type}_{mut}_LFC')
                if len(taa_naa_LFC_vals) == 0:
                    taa_wise_norm_LFC.append('-')
                else: 
                    taa_wise_norm_LFC.append(str(round(function_3Daggr(taa_naa_LFC_vals), 3)))
            df_struct_3d[f"{screen_name}_LFC3D"] = taa_wise_norm_LFC
        
        df_struct_3d[f"{screen_name}_LFC"] = df_struc_edits[f'{function_type}_{mut}_LFC']
        df_struct_3d[f"{screen_name}_LFC_Z"] = df_struc_edits[f'{function_type}_{mut}_LFC_Z']
    
        # GET RANDOMIZED LFC and LFC3D VALUES #
        if not os.path.exists(edits_filedir / rand_filename): 
            warnings.warn(f"{rand_filename} does not exist")
        df_struc_edits_rand = pd.read_csv(edits_filedir / rand_filename, sep = "\t")

        dict_temp = {}
        for r in range(0, nRandom):

            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = df_struc_edits_rand[f'{function_type}_missense_LFCr{str(r+1)}']

            if not LFC_only: 
                taa_wise_norm_LFC = []
                for aa in range(0, len(df_struc_edits_rand)):
                    taa_naa_LFC_vals = helper(df_struc_edits_rand, aa, lookup=f'{function_type}_missense_LFCr{str(r+1)}') ### issue this isn't randomized, is that ok?
                    if len(taa_naa_LFC_vals) == 0:
                        taa_wise_norm_LFC.append('-')
                    else:
                        taa_wise_norm_LFC.append(round(function_3Daggr(taa_naa_LFC_vals), 3))
                dict_temp[f"{screen_name}_LFC3Dr{str(r+1)}"] = taa_wise_norm_LFC

        df_struct_3d = pd.concat((df_struct_3d, pd.DataFrame(dict_temp)), axis=1)
        df_struct_3d = df_struct_3d.replace('-', np.nan)
        LFC_colnames   = [f"{screen_name}_LFCr{str(r+1)}" for r in range(0, nRandom)]
        LFC3D_colnames = [f"{screen_name}_LFC3Dr{str(r+1)}" for r in range(0, nRandom)]

        df_struct_3d = df_struct_3d.apply(lambda col: pd.to_numeric(col, errors='coerce'))
        df_struct_3d[f"{screen_name}_AVG_LFCr"]     = df_struct_3d[LFC_colnames].mean(axis=1) # AVG ALL
        df_struct_3d[f"{screen_name}_AVG_LFCr_neg"] = (df_struct_3d[LFC_colnames]
                                                        .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                        .sum(axis=1) / nRandom) # AVG NEG
        df_struct_3d[f"{screen_name}_AVG_LFCr_pos"] = (df_struct_3d[LFC_colnames]
                                                        .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                        .sum(axis=1) / nRandom) # AVG POS
        df_struct_3d = df_struct_3d.drop(LFC_colnames, axis=1)
        df_struct_3d[f"{screen_name}_AVG_LFC3Dr"]     = df_struct_3d[LFC3D_colnames].mean(axis=1) # AVG ALL
        df_struct_3d[f"{screen_name}_AVG_LFC3Dr_neg"] = (df_struct_3d[LFC3D_colnames]
                                                        .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                        .sum(axis=1) / nRandom) # AVG NEG
        df_struct_3d[f"{screen_name}_AVG_LFC3Dr_pos"] = (df_struct_3d[LFC3D_colnames]
                                                        .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                        .sum(axis=1) / nRandom) # AVG POS
        df_struct_3d = df_struct_3d.drop(LFC3D_colnames, axis=1)
        
        df_struct_3d = df_struct_3d.round(4)
        df_struct_3d = df_struct_3d.fillna('-')

    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr.tsv"
    df_struct_3d.to_csv(out_filename, sep = '\t', index=False)

    return df_struct_3d

def helper(
    df_struc_edits, aa, lookup
): 
    # naa IS NEIGHBORING AMINO ACIDS #
    # taa IS THIS AMINO ACID #
    taa_naa_LFC_vals = []
    taa_LFC = df_struc_edits.at[aa, lookup] # target LFC
    naa_pos_str = df_struc_edits.at[aa, 'Naa_pos']

    if taa_LFC != '-': # VALUE FOR THIS RESIDUE #
        taa_naa_LFC_vals.append(float(taa_LFC))

    if isinstance(naa_pos_str, str): # CHECK NEIGHBORING RESIDUES #
        naa_pos_list = naa_pos_str.split(';') # neighboring residue positions
        for naa_pos in naa_pos_list: 
            naa_LFC = df_struc_edits.at[int(naa_pos)-1, lookup]
            if naa_LFC != '-': 
                taa_naa_LFC_vals.append(float(naa_LFC))       

    return taa_naa_LFC_vals
