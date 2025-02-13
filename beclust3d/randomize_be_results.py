"""
File: randomized_preprocessed_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1.5

"""

import pandas as pd
from pathlib import Path
import numpy as np
import os

def randomize_be_results(df_missense, 
                         workdir, input_gene, screen_name, 
                         nRandom=1000, val_colname = 'LFC', 
                        #  seed=None, 
                         ): 
    """
    Description
        Takes reformatted missense dataframe and randomizes them to provide a baseline signal. 
        This function is run per gene per screen. 

    Params
        df_missense: pandas dataframe, required
            a dataframe listing all missense mutations and their LFC values
            with column headers edit, human_pos, refAA, altAA, LFC
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input gene ie 'TP53'
        input_screen: str, required
            the name of the input screen file ie 'ABE_screen_K562s.txt'
        screen_names: str, optional
            the name of the input screens
        nRandom: int, optional
            the number of randomize iterations
        val_colname: str, optional
            the name of the column with LFC values

    Returns
        df_missense: pandas dataframe
            a dataframe listing missense mutations with additional columns of randomized LFC values
            with column headers edit, human_pos, refAA, altAA, LFC, LFCr1 ... LFCr{nRandom}
    """
    # rng = np.random.default_rng(i)  # Fixed seed

    # NAME VARIABLES, PATHS, CREATE DIRECTORIES #
    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'randomized_screendata'):
        os.mkdir(edits_filedir / 'randomized_screendata')

    # SHUFFLE AND ADD TO DATAFRAME #
    LFC_list = df_missense[val_colname].tolist()
    dict_temp = {} # use a dictionary which is more efficient
    for i in range(nRandom): 
        dict_temp[f"{val_colname}r{str(i+1)}"] = np.random.permutation(LFC_list)
    df_missense = pd.concat((df_missense, pd.DataFrame(dict_temp)), axis=1)

    # SAVE RESULTS #
    out_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name.replace(' ','_')}_Missense_rand.tsv"
    df_missense.to_csv(out_filename, sep='\t', index=False)
    
    return df_missense
