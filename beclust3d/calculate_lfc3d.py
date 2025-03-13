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
        workdir, input_gene, screen_names, str_cons_dfs, str_cons_rand_dfs, 
        nRandom=1000, function_type='mean', mut='Missense', function_3Daggr=np.mean, 
        LFC_only=False, conserved_only=False, 
        # THERE ARE 2 MEAN FUNCTIONS, MEAN FOR CALCULATING LFC3D WHICH IS TUNABLE, AND MEAN FOR AVG RANDOMIZATIONS WHICH IS NOT TUNABLE #
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
    for screen_name, df_struc_edits, df_struc_edits_rand in zip(screen_names, str_cons_dfs, str_cons_rand_dfs):
        print(screen_name)

        if not LFC_only: # CALCULATE LFC3D #
            taa_wise_norm_LFC = []
            df_struc_edits_dict = df_struc_edits['conservation'].to_dict()
            taa_LFC_dict = df_struc_edits[f'{function_type}_{mut}_LFC'].to_dict()
            naa_pos_str_dict = df_struc_edits['Naa_pos'].to_dict()

            for aa in range(len(df_struc_edits)): # FOR EVERY RESIDUE #
                if conserved_only and df_struc_edits_dict[aa] != 'conserved': 
                    taa_wise_norm_LFC.append('-')
                    continue
                taa_naa_LFC_vals = helper(taa_LFC_dict, df_struc_edits_dict, aa, 
                                          naa_pos_str_dict[aa], conserved_only)
                if len(taa_naa_LFC_vals) == 0:
                    taa_wise_norm_LFC.append('-')
                else: 
                    taa_wise_norm_LFC.append(str(function_3Daggr(taa_naa_LFC_vals)))
            df_struct_3d = pd.concat([df_struct_3d, pd.DataFrame({f"{screen_name}_LFC3D": taa_wise_norm_LFC})], axis=1)
            del df_struc_edits_dict, taa_wise_norm_LFC
        
        # df_struct_3d[f"{screen_name}_LFC"] = df_struc_edits[f'{function_type}_{mut}_LFC'] # duplicate LFC
        # df_struct_3d[f"{screen_name}_LFC_Z"] = df_struc_edits[f'{function_type}_{mut}_LFC_Z'] # duplicate LFC-Z

        df_struct_3d = pd.concat([df_struct_3d, 
                                  df_struc_edits[[f'{function_type}_{mut}_LFC']].rename(
                                      columns={f'{function_type}_{mut}_LFC': f"{screen_name}_LFC"}), 
                                  df_struc_edits[[f'{function_type}_{mut}_LFC_Z']].rename(
                                      columns={f'{function_type}_{mut}_LFC_Z': f"{screen_name}_LFC_Z"}), 
                                  ], axis=1)

        dict_temp = {}
        for r in range(0, nRandom):

            dict_temp[f"{screen_name}_LFCr{str(r+1)}"] = df_struc_edits_rand[f'{function_type}_missense_LFCr{str(r+1)}']

            if not LFC_only: 
                taa_wise_norm_LFC = []
                df_struc_edits_rand_dict = df_struc_edits_rand['conservation'].to_dict()
                taa_LFC_dict = df_struc_edits_rand[f'{function_type}_missense_LFCr{str(r+1)}'].to_dict()
                naa_pos_str_dict = df_struc_edits_rand['Naa_pos'].to_dict()

                for aa in range(len(df_struc_edits_rand)):
                    if conserved_only and df_struc_edits_rand_dict[aa] != 'conserved': 
                        taa_wise_norm_LFC.append('-')
                        continue
                    # taa_naa_LFC_vals = helper(df_struc_edits_rand, aa, f'{function_type}_missense_LFCr{str(r+1)}', conserved_only) 
                    taa_naa_LFC_vals = helper(taa_LFC_dict, df_struc_edits_rand_dict, aa, 
                                              naa_pos_str_dict[aa], conserved_only)
                    if len(taa_naa_LFC_vals) == 0:
                        taa_wise_norm_LFC.append('-')
                    else:
                        taa_wise_norm_LFC.append(function_3Daggr(taa_naa_LFC_vals))
                dict_temp[f"{screen_name}_LFC3Dr{str(r+1)}"] = taa_wise_norm_LFC
                del df_struc_edits_rand_dict, taa_wise_norm_LFC

        df_struct_3d = pd.concat((df_struct_3d, pd.DataFrame(dict_temp)), axis=1)
        df_struct_3d = df_struct_3d.replace('-', np.nan).infer_objects(copy=False) # FutureWarning
        LFC_colnames   = [f"{screen_name}_LFCr{str(r+1)}" for r in range(0, nRandom)]
        del dict_temp

        df_struct_3d = df_struct_3d.apply(lambda col: pd.to_numeric(col, errors='coerce'))
        df_struct_3d[f"{screen_name}_AVG_LFCr"]     = df_struct_3d[LFC_colnames].mean(axis=1) # AVG ALL
        df_struct_3d[f"{screen_name}_AVG_LFCr_neg"] = (df_struct_3d[LFC_colnames]
                                                        .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                        .sum(axis=1) / nRandom) # AVG NEG
        df_struct_3d[f"{screen_name}_AVG_LFCr_pos"] = (df_struct_3d[LFC_colnames]
                                                        .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                        .sum(axis=1) / nRandom) # AVG POS
        # df_struct_3d = df_struct_3d.drop(LFC_colnames, axis=1)
        if not LFC_only: 
            LFC3D_colnames = [f"{screen_name}_LFC3Dr{str(r+1)}" for r in range(0, nRandom)]
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr"]     = df_struct_3d[LFC3D_colnames].mean(axis=1) # AVG ALL
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr_neg"] = (df_struct_3d[LFC3D_colnames]
                                                            .apply(lambda col: col.map(lambda x: x if x < 0 else np.nan))
                                                            .sum(axis=1) / nRandom) # AVG NEG
            df_struct_3d[f"{screen_name}_AVG_LFC3Dr_pos"] = (df_struct_3d[LFC3D_colnames]
                                                            .apply(lambda col: col.map(lambda x: x if x > 0 else np.nan))
                                                            .sum(axis=1) / nRandom) # AVG POS
            # df_struct_3d = df_struct_3d.drop(LFC3D_colnames, axis=1)
        
        df_struct_3d = df_struct_3d
        df_struct_3d = df_struct_3d.fillna('-')

    df_struct_3d['unires'] = df_str_cons['unires']
    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_LFC3Dr.tsv"
    df_struct_3d.to_csv(out_filename, sep = '\t', index=False)

    return df_struct_3d

def helper(
    taa_LFC_dict, df_struc_edits_dict, aa, 
    naa_pos_str, conserved_only
): 
    # naa IS NEIGHBORING AMINO ACIDS #
    # taa IS THIS AMINO ACID #
    taa_naa_LFC_vals = []
    taa_LFC = taa_LFC_dict[aa]

    if taa_LFC != '-': # VALUE FOR THIS RESIDUE #
        if not conserved_only or df_struc_edits_dict[aa] == 'conserved': 
            taa_naa_LFC_vals.append(float(taa_LFC))

    if isinstance(naa_pos_str, str): # CHECK NEIGHBORING RESIDUES #
        naa_pos_list = naa_pos_str.split(';') # neighboring residue positions
        for naa_pos in naa_pos_list: 
            if not conserved_only or df_struc_edits_dict[int(naa_pos)-1] == 'conserved': 
                naa_LFC = taa_LFC_dict[int(naa_pos)-1]
                if naa_LFC != '-': 
                    taa_naa_LFC_vals.append(float(naa_LFC))

    return taa_naa_LFC_vals
