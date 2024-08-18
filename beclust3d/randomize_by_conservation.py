"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Translated from Notebook 3.2.5

"""

import pandas as pd
from pathlib import Path
import os
import warnings
import numpy as np

def randomize_by_conservation(
        workdir, 
        input_gene, input_screen, structureid, 
        nRandom=1000, 
): 
    """
    Description
        Randomizes the scores weighted by structural conservation fom previous step

    Params
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        structureid: str, required
            the name of the AF and uniprot input
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_struc_consvr: pandas dataframe
            a dataframe listing randomized structural conservation data
    """

    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir)
    edits_filedir = edits_filedir / input_gene
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'randomized_screendata'):
        os.mkdir(edits_filedir / 'randomized_screendata')
    
    struc_consrv_filename =  f"screendata/{input_gene}_{structureid}_struc_consrv.tsv"
    if not os.path.exists(edits_filedir / struc_consrv_filename): 
        warnings.warn(f"{missense_filename} does not exist")
        return None
    df_struc_consvr = pd.read_csv(edits_filedir / struc_consrv_filename, sep = "\t")

    missense_filename = f"randomized_screendata/{input_gene}_{screen_name}_missense_edits_randomized.tsv"
    if not os.path.exists(edits_filedir / missense_filename): 
        warnings.warn(f"{missense_filename} does not exist")
        return None
    df_missense = pd.read_csv(edits_filedir / missense_filename, sep = '\t') # DOESNT CONTAIN '-' #

    # COPY VALUES FROM MISSENSE RANDOMIZED DF TO STRUC CONSRV DF #
    human_res_positions = df_struc_consvr['human_res_pos'].tolist()
    missense_filter_col = [col for col in df_missense if col.startswith('LFC')]
    
    new_colnames = ['mean_missense_LFC']+[f'mean_missense_LFCr{j+1}' for j in range(nRandom)] # +['mean_missense_LFCavg']
    df_mis_positions = pd.DataFrame(columns=new_colnames)
    for i in range(len(df_struc_consvr)): 
        df_mis_pos = df_missense.loc[df_missense['edit_pos'] == human_res_positions[i]]
        df_mis_pos = df_mis_pos[missense_filter_col]

        if df_mis_pos.shape[0] == 0: # FILL WITH '-' #
            res = ['-' for _ in range(df_mis_pos.shape[1])]
        else: # AVERAGE ACROSS COLUMNS FOR ONE OR MORE ROWS #
            res = df_mis_pos.mean().tolist()
        df_mis_positions.loc[i] = res # ADD ROW FROM MISSENSE RANDOMIZED #

    df_mis_positions = df_mis_positions.round(3)
    df_struc_consvr = pd.concat([df_struc_consvr, df_mis_positions], axis=1)
    del df_mis_pos, df_mis_positions, df_missense

    out_filename = f"randomized_screendata/{input_gene}_{screen_name}_struc_consrv_missenseedits_randomized.tsv"
    df_struc_consvr.to_csv(edits_filedir / out_filename, sep = '\t', index=False)

    return df_struc_consvr
