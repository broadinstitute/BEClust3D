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
warnings.filterwarnings('ignore')

def randomize_by_sequence(
        workdir, 
        df_missense, df_rand, 
        input_gene, screen_name, 
        nRandom=1000, conservation=False, 
        function_name='mean', target_pos='human_res_pos', target_res='', 
        # THERE ARE 2 MEAN FUNCTIONS, MEAN FOR CALCULATING LFC3D WHICH IS TUNABLE, AND MEAN FOR AVG RANDOMIZATIONS WHICH IS NOT TUNABLE #
): 
    """
    Description
        Randomizes the scores weighted by structural sequence and conservation fom previous step

    Params
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_screen: str, required
            the name of the input screen
        screen_names: str, optional
            the name of the input screen
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_missense: pandas dataframe
            a dataframe listing randomized structural sequence and conservation data
    """

    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'randomized_screendata'):
        os.mkdir(edits_filedir / 'randomized_screendata')

    # COPY VALUES FROM MISSENSE RANDOMIZED DF TO STRUC CONSRV DF #
    res_position_list = df_missense[target_pos].tolist()
    if len(target_res) > 0 and target_res in df_missense.columns: # GET RID OF NON CONSERVED FOR MOUSE RES ONLY #
        res_list = df_missense[target_res].tolist()
    else: res_list = None
    missense_filter_col = [col for col in df_rand.columns if col.startswith('LFC')]
    
    rand_colnames = [f'{function_name}_missense_LFC']+[f'{function_name}_missense_LFCr{j+1}' for j in range(nRandom)]
    df_mis_positions = pd.DataFrame(columns=rand_colnames)
    for i in range(len(df_missense)): 
        df_mis_pos = df_rand.loc[df_rand['edit_pos'] == res_position_list[i]]
        df_mis_pos = df_mis_pos[missense_filter_col]

        if df_mis_pos.shape[0] == 0 and (res_list is not None and res_list[i] != '-'): # FILL WITH '-' #
            res = ['-' for _ in range(df_mis_pos.shape[1])]
        else: # AVERAGE ACROSS COLUMNS FOR ONE OR MORE ROWS #
            res = df_mis_pos.mean().tolist()
        df_mis_positions.loc[i] = res # ADD ROW FROM MISSENSE RANDOMIZED #

    missense_colnames = ['unipos', 'unires', 'x_coord', 'y_coord', 'z_coord', 
                         'bfactor_pLDDT', 'Naa_count', 'Naa', 'Naa_pos', 'SS9', 'SS3', 'ACC', 
                         'RSA', 'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI', 'dBurial', 
                         'normSumdBurial', 'pLDDT_dis', 'human_res_pos', 'conservation', 
                         f'{function_name}_Missense_LFC',f'{function_name}_Missense_LFC_stdev','all_Missense_edits',
                         f'{function_name}_Missense_LFC_Z',f'{function_name}_Missense_LFC_p',f'{function_name}_Missense_LFC_plab']
    if conservation: 
        missense_colnames += ['mouse_res_pos', 'mouse_res']
    df_mis_positions = df_mis_positions
    df_missense = pd.concat([df_missense[missense_colnames], df_mis_positions], axis=1)
    df_missense = df_missense
    del df_mis_positions, df_rand

    out_filename = f"randomized_screendata/{input_gene}_{screen_name}_Missense_proteinedits_rand.tsv"
    df_missense.to_csv(edits_filedir / out_filename, sep = '\t', index=False)

    return df_missense
