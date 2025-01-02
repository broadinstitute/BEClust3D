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

def randomize_by_sequence(
        workdir, 
        struc_consrv_filename, missense_filename, 
        input_gene, screen_name, 
        nRandom=1000, 
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
        df_protein: pandas dataframe
            a dataframe listing randomized structural sequence and conservation data
    """

    edits_filedir = Path(workdir)
    edits_filedir = edits_filedir / input_gene
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'randomized_screendata'):
        os.mkdir(edits_filedir / 'randomized_screendata')
    
    if not os.path.exists(edits_filedir / struc_consrv_filename): 
        warnings.warn(f"{struc_consrv_filename} does not exist")
        return None
    df_protein = pd.read_csv(edits_filedir / struc_consrv_filename, sep = "\t")

    if not os.path.exists(edits_filedir / missense_filename): 
        warnings.warn(f"{missense_filename} does not exist")
        return None
    df_missense = pd.read_csv(edits_filedir / missense_filename, sep = '\t') # DOESNT CONTAIN '-' #

    # COPY VALUES FROM MISSENSE RANDOMIZED DF TO STRUC CONSRV DF #
    human_res_positions = df_protein['human_res_pos'].tolist()
    missense_filter_col = [col for col in df_missense if col.startswith('LFC')]
    
    new_colnames = ['mean_missense_LFC']+[f'mean_missense_LFCr{j+1}' for j in range(nRandom)]
    df_mis_positions = pd.DataFrame(columns=new_colnames)
    for i in range(len(df_protein)): 
        df_mis_pos = df_missense.loc[df_missense['edit_pos'] == human_res_positions[i]]
        df_mis_pos = df_mis_pos[missense_filter_col]

        if df_mis_pos.shape[0] == 0: # FILL WITH '-' #
            res = ['-' for _ in range(df_mis_pos.shape[1])]
        else: # AVERAGE ACROSS COLUMNS FOR ONE OR MORE ROWS #
            res = df_mis_pos.mean().tolist()
        df_mis_positions.loc[i] = res # ADD ROW FROM MISSENSE RANDOMIZED #

    df_mis_positions = df_mis_positions.round(3)
    df_protein = pd.concat([df_protein, df_mis_positions], axis=1)
    del df_mis_pos, df_mis_positions, df_missense

    out_filename = f"randomized_screendata/{input_gene}_{screen_name}_Missense_proteinedits_rand.tsv"
    df_protein.to_csv(edits_filedir / out_filename, sep = '\t', index=False)

    return df_protein
