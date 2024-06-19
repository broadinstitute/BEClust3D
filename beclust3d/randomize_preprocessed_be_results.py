"""
File: randomized_preprocessed_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1.5

"""

import pandas as pd
from pathlib import Path
import random 

def randomize_be_results(df_missense, workdir, 
                         input_gene, input_screen, 
                         nRandom=1000
                         ): 
    """
    Description
        Takes reformatted missense dataframe and randomizes them to provide a baseline signal. 

    Params

    Returns
    
    """

    ### figure out creating directories

    edits_filedir = Path(workdir + input_gene)

    LFC_list = df_missense["LFC"].tolist()
    for i in range(0, nRandom):
        LFC_list_random = random.sample(LFC_list, len(LFC_list))
        headertext_new = "LFCr" + str(i+1)
        df_missense[headertext_new] = LFC_list_random

    # save results
    out_filename = edits_filedir / f"randomized_screendata/{input_gene}_{input_screen}_missense_edits_randomized.tsv"
    df_missense.to_csv(out_filename, sep='\t', index=False)
    
    return df_missense

