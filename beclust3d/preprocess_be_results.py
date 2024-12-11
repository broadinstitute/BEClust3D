"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1

"""

import re
from pathlib import Path
import numpy as np
import os
import re
import warnings

from _preprocess_be_results_plots_ import *

aa_map = {
    'ALA': 'A',  'ARG': 'R',  'ASN': 'N',  'ASP': 'D',  'CYS': 'C',  
    'GLN': 'Q',  'GLU': 'E',  'GLY': 'G',  'HIS': 'H',  'ILE': 'I',  
    'LEU': 'L',  'LYS': 'K',  'MET': 'M',  'PHE': 'F',  'PRO': 'P',  
    'SER': 'S',  'THR': 'T',  'TRP': 'W',  'TYR': 'Y',  'VAL': 'V', 
}

# These have to be predefined, too much of the code is dependent on these categories #
mut_categories = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"]

# change df_Input and input_screen in lists
# iterate through both and aggregate

def parse_base_editing_results(
    input_dfs, workdir, 
    input_gene, input_screens, screen_names=[], 
    mut_col='Mutation category', val_col='logFC', 
    gene_col='Target Gene Symbol', edits_col='Amino Acid Edits', 
    multi_annotation = False, split_char = ',',
): 
    """
    Description
        Parse raw data and create separate dataframes for each mutation type. 

    Params
        df_Input: pandas dataframe
            the raw dataset containing columns mut_col, val_col, gene_col, edits_col
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input gene
        input_screen: list of str, required
            the names of the input screen files
        screen_names: list of str, optional
            the names of the input screens
        mut_col: str, optional
            column name that indicates the type of mutation
        val_col: str, optional
            column name that indicates the log FC value
        gene_col: str, optional
            column name that indicates the name of the gene
        edits_col: str, optional
            column name that indicates the list of edits
        multi_annotation: bool, optional
            whether or not mut_col is list of mutations rather than one mutation
        split_char: str, optional
            the delimiter on which to split the edits (ie T78A,F79P)

    Returns
        mut_dfs: a dictionary of pandas dataframes
            each dataframe corresponds to data for Missense, Silent, Nonsense, Intron, UTR, No mutation
            each dataframe with column headers edit, {output_col}, refAA, altAA, LFC
    """

    edits_filedir = Path(workdir)
    edits_filedir = edits_filedir / input_gene
    if not screen_names: # for screen_names being an empty list
        screen_names = [input_screen.split('.')[0] for input_screen in input_screens]
    if not os.path.exists(edits_filedir): 
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')
    if not os.path.exists(edits_filedir / 'qc_validation'):
        os.mkdir(edits_filedir / 'qc_validation')
    
    # INDIVIDUAL BARPLOTS AND VIOLIN PLOTS FOR EACH SCREEN #
    for df, screen_name in zip(input_dfs, screen_names): 
        counts_by_gene(df_inputs=[df], edits_filedir=edits_filedir, 
                       gene_col=gene_col, mut_col=mut_col, title=screen_name)
        violin_by_gene(df_inputs=[df], edits_filedir=edits_filedir, 
                       gene_col=gene_col, mut_col=mut_col, val_col=val_col, title=screen_name)
    # AGGREGATE ACROSS SCREENS FOR SUMMARY PLOTS #
    if len(input_dfs) > 1: 
        counts_by_gene(df_inputs=input_dfs, edits_filedir=edits_filedir, 
                       gene_col=gene_col, mut_col=mut_col, title='Aggregate')
        violin_by_gene(df_inputs=input_dfs, edits_filedir=edits_filedir, 
                       gene_col=gene_col, mut_col=mut_col, val_col=val_col, title='Aggregate')

    mut_dfs = {}
    # OUTPUT TSV BY INDIVIDUAL SCREENS #
    for df_Input, screen_name in zip(input_dfs, screen_names): 
        # NARROW DOWN TO INPUT_GENE #
        df_InputGene = df_Input.loc[df_Input[gene_col] == input_gene, ]
        mut_dfs[screen_name] = {}
        for mut_cat in mut_categories: 
            if not mut_cat in df_InputGene[mut_col].unique(): 
                warnings.warn(f'{mut_cat} not in Dataframe')

        # NARROW DOWN TO EACH MUTATION TYPE #
        for mut in mut_categories: 

            # IF USER WANTS TO CATEGORIZE BY ONE SINGLE MUTATION PER GUIDE OR MULTIPLE MUTATIONS PER GUIDE #
            if multi_annotation: 
                df = df_InputGene.loc[df_InputGene[mut_col].contains(mut), ]
            else: 
                df = df_InputGene.loc[df_InputGene[mut_col] == mut, ]
            df = df.reset_index(drop=True)
            print(f"Count of {mut} rows: " + str(len(df)))

            # ASSIGN position refAA altAA #
            df[edits_col] = df[edits_col].str.strip(',').str.strip(';') # CLEAN
            df[edits_col] = df[edits_col].str.split(split_char) # STR to LIST
            df[val_col] = df[val_col].round(3)
            df[edits_col] = df[edits_col].apply(lambda xs: [x for x in xs if re.match('^[A-Z*][0-9]{1,4}[A-Z*]$', x)] if not isinstance(xs, float) else []) # FILTER FOR MUTATIONS #

            df_exploded = df.explode(edits_col) # EACH ROW IS A MUTATION #
            df_exploded['edit_pos'] = df_exploded[edits_col].str.extract('(\d+)')
            df_exploded['refAA'], df_exploded['altAA'] = df_exploded[edits_col].str[0], df_exploded[edits_col].str[-1]
            # IF 3 LETTER CODES ARE USED, TRANSLATE TO 1 LETTER CODE #
            df_exploded['refAA'] = df_exploded['refAA'].str.upper().map(aa_map)
            df_exploded['altAA'] = df_exploded['altAA'].str.upper().map(aa_map)
            df_subset = df_exploded[[edits_col, 'edit_pos', 'refAA', 'altAA', val_col]].rename(columns={edits_col: 'this_edit', val_col: 'LFC'})

            if mut == 'Missense': 
                df_subset = df_subset[(df_subset['refAA'] != df_subset['altAA']) & (df_subset['altAA'] != '*')]
            elif mut == 'Silent': # SILENT BEFORE NONSENSE (ie *248* MUTATION IS SILENT NOT NONSENSE)
                df_subset = df_subset[df_subset['refAA'] == df_subset['altAA']]
            elif mut == 'Nonsense': 
                df_subset = df_subset[df_subset['altAA'] == '*']
            else: 
                df_subset = df_subset['LFC']

            # WRITE LIST OF MUT AND THEIR LFC VALUES #
            edits_filename = f"screendata/{input_gene}_{screen_name.replace(' ','_')}_{mut.replace(' ','_')}.tsv"
            df_subset.to_csv(edits_filedir / edits_filename, sep='\t')
            mut_dfs[screen_name][mut] = df_subset
    
    # AGGREGATE ACROSS SCREENS FOR PLOTS #
    # MANN WHITNEY TEST and VIOLIN PLOTS #
    if len(input_screens) > 1:
        for df, screen_name in zip(input_dfs, screen_names): 
            df_muts, mw_res = mann_whitney_test(edits_filedir, [screen_name], input_gene)
            df_muts['LFC_direction'] = np.where(df_muts['LFC'] < 0, 'neg', 'pos')
            # VIOLIN PLOTS #
            violin_plot(df_muts, edits_filedir, input_gene, screen_name)
        if len(input_dfs) > 1: 
            df_muts, mw_res = mann_whitney_test(edits_filedir, screen_names, input_gene, )
            df_muts['LFC_direction'] = np.where(df_muts['LFC'] < 0, 'neg', 'pos')
            # VIOLIN PLOTS #
            violin_plot(df_muts, edits_filedir, input_gene, 'Aggregate')

    return mut_dfs
