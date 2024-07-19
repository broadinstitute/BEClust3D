"""
File: randomized_preprocessed_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1.5

"""

import pandas as pd
from pathlib import Path
import random
import os

def randomize_be_results(df_missense, workdir, 
                         input_gene, input_screen, 
                         nRandom=1000
                         ): 
    """
    Description
        Takes reformatted missense dataframe and randomizes them to provide a baseline signal. 

    Params
        df_missense: pandas dataframe, required
            a dataframe listing all missense mutations and their LFC values
            with column headers edit, human_pos, refAA, altAA, LFC
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input gene
        input_screen: str, required
            the name of the input screen
        nRandom: int, optional
            the number of randomize iterations

    Returns
        df_missense: pandas dataframe
            a dataframe listing missense mutations with additional columns of randomized LFC values
            with column headers edit, human_pos, refAA, altAA, LFC, LFCr1 ... LFCr{nRandom}
    """

    screen_name = input_screen.split('.')[0]
    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'randomized_screendata'):
        os.mkdir(edits_filedir / 'randomized_screendata')

    LFC_list = df_missense["LFC"].tolist()
    dict_temp = {}
    for i in range(0, nRandom):
        LFC_list_random = random.sample(LFC_list, len(LFC_list))
        headertext_new = "LFCr" + str(i+1)
        dict_temp[headertext_new] = LFC_list_random
    # Using pd.concat does not throw a performance warning, but takes 6x the time if I pd.concat every step
    df_missense = pd.concat((df_missense, pd.DataFrame(dict_temp)), axis=1) # Use dict instead

    # save results
    out_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screen_name}_missense_edits_randomized.tsv"
    df_missense.to_csv(out_filename, sep='\t', index=False)
    
    return df_missense
