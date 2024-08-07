"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1

"""

import pandas as pd
import re
from pathlib import Path
import numpy as np
import os
import re

from _preprocess_be_results_plots_ import *
from _preprocess_be_results_hypothesis_ import *

# These have to be predefined, too much of the code is dependent on these categories #
mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"]

def parse_base_editing_results(
    df_Input, workdir, 
    input_gene, input_screen, 
    mut_col='Mutation category', val_col='logFC', 
    gene_col='Target Gene Symbol', edits_col='Amino Acid Edits', 
    output_col='human_pos', 
    multi_annotation = False, split_char = ',',
): 
    """
    Description
        Parse raw data and create separate dataframes for each mutation type. 

    Params
        df_InputGene: pandas dataframe
            the raw dataset containing columns mut_col, val_col, gene_col, edits_col
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input gene
        input_screen: str, required
            the name of the input screen
        mut_col: str, optional
            column name that indicates the type of mutation
        val_col: str, optional
            column name that indicates the log FC value
        gene_col: str, optional
            column name that indicates the name of the gene
        edits_col: str, optional
            column name that indicates the list of edits

    Returns
        mut_dfs: a list of pandas dataframes
            each dataframe corresponds to data for Missense, Silent, Nonsense, Intron, UTR, No mutation
            each dataframe with column headers edit, {output_col}, refAA, altAA, LFC
    """

    edits_filedir = Path(workdir + '/' +  input_gene)
    screen_name = input_screen.split('.')[0]    
    counts_violin_by_gene(df_rawinput=df_Input, edits_filedir=edits_filedir, 
                          gene_col=gene_col, mut_col=mut_col, val_col=val_col, screen_name=screen_name, )
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')
    if not os.path.exists(edits_filedir / 'qc_validation'):
        os.mkdir(edits_filedir / 'qc_validation')

    # NARROW DOWN TO INPUT_GENE #
    df_InputGene = df_Input.loc[df_Input[gene_col] == input_gene, ]
    mut_dfs = []
    for mut_cat in mut_categories: 
        assert mut_cat in df_InputGene[mut_col].unique()

    # NARROW DOWN TO EACH MUTATION TYPE #
    for mut in mut_categories: 
        df = pd.DataFrame()
        # IF USER WANTS TO CATEGORIZE BY ONE SINGLE MUTATION PER GUIDE OR MULTIPLE MUTATIONS PER GUIDE #
        if multi_annotation: 
            df = df_InputGene.loc[df_InputGene[mut_col].contains(mut), ]
        else: 
            df = df_InputGene.loc[df_InputGene[mut_col] == mut, ]
        df = df.reset_index(drop=True)
        print(f"Count of {mut} rows: " + str(len(df)))

        # ASSIGN position refAA altAA #
        df[edits_col] = df[edits_col].str.strip(',').str.strip(';').str.split(split_char)
        df[val_col] = df[val_col].round(3)
        df[edits_col] = df[edits_col].apply(lambda xs: [x for x in xs if re.match('^[A-Z*][0-9]{1,4}[A-Z*]$', x)] if not isinstance(xs, float) else []) # FILTER FOR MUTATIONS #

        df_exploded = df.explode(edits_col) # EACH ROW IS A MUTATION #
        df_exploded['edit_pos'] = df_exploded[edits_col].str.extract('(\d+)')
        df_exploded['refAA'], df_exploded['altAA'] = df_exploded[edits_col].str[0], df_exploded[edits_col].str[-1]
        df_subset = df_exploded[[edits_col, 'edit_pos', 'refAA', 'altAA', val_col]].rename(columns={edits_col: 'this_edit', val_col: 'LFC'})

        if mut == 'Missense': 
            df_subset = df_subset[(df_subset['refAA'] != df_subset['altAA']) & (df_subset['altAA'] != '*')]
        elif mut == 'Nonsense': 
            df_subset = df_subset[df_subset['altAA'] == '*']
        elif mut == 'Silent': 
            df_subset = df_subset[df_subset['refAA'] == df_subset['altAA']]
        else: 
            df_subset = df_subset['LFC']

        # WRITE LIST OF MUT AND THEIR LFC VALUES #
        edits_filename = f"screendata/{input_gene}_{screen_name}_{mut.replace(' ', '_')}_edits_list.tsv"
        df_subset.to_csv(edits_filedir / edits_filename, sep='\t')
        mut_dfs.append(df_subset)
    
    # MW AND KS TESTS HYPOTHESIS 1 #
    hypothesis_one(df_Input, edits_filedir, 
                   screen_name, gene_col, mut_col, val_col, testtype='MannWhitney')
    hypothesis_one(df_Input, edits_filedir, 
                   screen_name, gene_col, mut_col, val_col, testtype='KolmogorovSmirnov')
    hypothesis_plot(edits_filedir, screen_name, testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='1')
    # MW AND KS TESTS HYPOTHESIS 2 #
    hypothesis_two(df_Input, edits_filedir, 
                   screen_name, gene_col, mut_col, val_col, testtype='MannWhitney')
    hypothesis_two(df_Input, edits_filedir, 
                   screen_name, gene_col, mut_col, val_col, testtype='KolmogorovSmirnov')
    hypothesis_plot(edits_filedir, screen_name, testtype1='MannWhitney', testtype2='KolmogorovSmirnov', hypothesis='2')
    
    # TESTS AND PLOTS #
    # MANN WHITNEY TEST #
    df_list, mw_res = mann_whitney_test(edits_filedir=edits_filedir, 
                                        screen_name=screen_name, input_gene=input_gene, )
    print(df_list)
    # VIOLIN PLOTS #
    violin_plot(df_InputGene_edits_list=df_list, 
                edits_filedir=edits_filedir, screen_name=screen_name, input_gene=input_gene, )
    # DIRECTIONAL VIOLIN PLOTS #
    df_list['LFC_direction'] = np.where(df_list['LFC'] < 0, 'neg', 'pos')
    violin_plot(df_InputGene_edits_list=df_list, 
                edits_filedir=edits_filedir, screen_name=screen_name, input_gene=input_gene, 
                directional=True, )

    return mut_dfs
