"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1

"""

import pandas as pd
import re
from pathlib import Path
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def parse_base_editing_results(df_InputGene, workdir, 
                               input_gene, input_screen, 
                               mut_col='Mutation category', val_col='logFC', 
                               gene_col='Target Gene Symbol', edits_col='Amino Acid Edits', 
                               output_col='human_pos', 
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

    ### figure out creating directories

    df_InputGene = df_InputGene.loc[df_InputGene[gene_col] == input_gene, ]
    screen_name = input_screen.split('.')[0]
    # df_InputGene_missense, df_InputGene_silent, df_InputGene_nonsense, df_InputGene_intron, df_InputGene_UTR
    mut_categories = ['Missense', 'Silent', 'Nonsense'] # 'Intron', 'UTR'
        # based on the standard output from BEAGLE
        # order of priority
    mut_dfs = [pd.DataFrame() for _ in mut_categories]

    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'screendata'):
        os.mkdir(edits_filedir / 'screendata')
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')

    for mut, df in zip(mut_categories, mut_dfs): 
        # preprocess by taking the subset of the original df with the indicated mutation
        df = df_InputGene.loc[df_InputGene[mut_col] == mut, ] ### == or .contains
        df = df.reset_index(drop=True)
        print(f"Count of {mut} rows: " + str(len(df)))

        # only Missense, Silent, Nonsense have associated mutations to parse
        if mut not in ['Missense', 'Silent', 'Nonsense']: 
            continue

        # open a file to write in list of mutations and their LFC values
        edits_filename = edits_filedir / f"screendata/{input_gene}_{screen_name}_{mut}_edits_list.tsv"
        edits_file = open(edits_filename, "w")
        edits_file.write("edit\t"+output_col+"\trefAA\taltAA\tLFC\n")

        # iterate through each cell of each list of edits
        for i in range(0, len(df)):
            edits_all = df.at[i, edits_col].strip(',').strip(';')
            edit_val = round(df.at[i, val_col], 3)
            edits_list = edits_all.split(',')
            
            for j in range(0, len(edits_list)):
                # identify edit information and output them into .tsv
                this_edit = edits_list[j].strip()
                temp = parse_edit_helper(mut, this_edit)
                if temp is not None: 
                    edit_refAA, edit_altAA, edit_pos = temp
                    line = '\t'.join([this_edit, str(edit_pos), edit_refAA, edit_altAA, str(edit_val)]) + '\n'
                    edits_file.write(line)
                
        edits_file.close()

    # preprocess these separately
    for mut in ['Splice Site', 'No Mutation']: 
        df_temp = df_InputGene.loc[df_InputGene[mut_col] == mut, ]
        df_temp = df_temp.reset_index(drop=True)
        print("Count of Splice site rows: " + str(len(df_temp)))

        df_temp_new = pd.DataFrame()
        df_temp_new['gene'] = df_temp[gene_col]
        df_temp_new['LFC'] = df_temp[val_col]

        mut_str = mut.replace(' ', '_')
        df_temp_new.to_csv(edits_filedir / f"screendata/{input_gene}_{screen_name}_{mut_str}_edits_list.tsv", 
                           sep='\t', index=False)
        mut_dfs.append(df_temp)

    # TESTS and PLOTS #
    # mann whitney test
    df_list, comparisons = mann_whitney_test(edits_filedir=edits_filedir, 
                                screen_name=screen_name, input_gene=input_gene)

    # violin plot by mutation
    violin_plot(
            df_InputGene_edits_list=df_list, 
            edits_filedir=edits_filedir, screen_name=screen_name, input_gene=input_gene,
            )
    # directional violin plot by mutation
    df_list['LFC_direction'] = np.where(df_list['LFC'] < 0, 'neg', 'pos')
    violin_plot(
            df_InputGene_edits_list=df_list, 
            edits_filedir=edits_filedir, screen_name=screen_name, input_gene=input_gene, 
            directional=True,
            )

    return mut_dfs


def parse_edit_helper(
        mut_type, this_edit
): 
    """
    Description
        A helper function to take in a string and parse out the mutation information. 
        For example, Met20Ala is interpreted as 'Met', 20, 'Ala'

    Params
        mut_type: str, required
            one of the types of mutations ie Missense, Silent, Nonsense, etc
        this_edit: str, required
            a string in the approximate format of Met20Ala or M20A

    Returns
        edit_refAA: str
            the original amino acid
        edit_altAA: str
            the new amino acid
        edit_pos: int
            the position of the base edit
    """
    ### different format types
    pattern = r'^([a-zA-Z*]{1,3})(\d{1,4})([a-zA-Z*]{1,3})$'
    match_edit = re.match(pattern, this_edit)
    if match_edit:
        edit_refAA, edit_pos, edit_altAA = list(match_edit.groups())

        ### what do i do if assertions fail
        if mut_type == "Missense" and edit_refAA != edit_altAA: 
            return edit_refAA, edit_altAA, edit_pos

        elif mut_type == "Silent" and edit_refAA == edit_altAA:
            return edit_refAA, edit_altAA, edit_pos

        elif mut_type == "Nonsense" and edit_altAA in ['Ter', 'STOP', '*']: ### some possible ways to decribe stop codon
            return edit_refAA, edit_altAA, edit_pos
    
    return None

def mann_whitney_test(
    edits_filedir, screen_name, input_gene
): 
    """
    Description
        A helper function to run the Mann Whitney test on
        'Missense', 'Silent', 'Nonsense', 'No_mutation'

    Params
        edits_filedir: Path, required
            Path to working directory
        screen_name: str, required
            the name of the input screen parsed form the input file
        input_gene: str, required
            the name of the input gene

    Returns
        df_InputGene_edits_list: list of Pandas Dataframes
            A list of dataframes, for each category of mutations
        comparisons: dict
            A dictionary of each screen comparison and their Mann Whitney results
    """
    
    mut_categories = ['Missense', 'Silent', 'Nonsense', 'No_mutation']
    df_inputgenes = [pd.DataFrame() for _ in mut_categories]

    for mut, df in zip(mut_categories, df_inputgenes): 
        filename = edits_filedir / f"screendata/{input_gene}_{screen_name}_{mut}_edits_list.tsv"
        df_temp = pd.read_csv(filename, sep = '\t')
        df['LFC'] = df_temp['LFC']
        df['muttype'] = mut
        print(f"{mut} edits: {str(len(df))}")

    # Mann Whitney U test #
    comparisons = [
        {'df1': df_inputgenes[2], 'df2': df_inputgenes[0], 'name': 'Nonsense vs Missense'}, 
        {'df1': df_inputgenes[2], 'df2': df_inputgenes[3], 'name': 'Nonsense vs No mutation'}, 
        {'df1': df_inputgenes[2], 'df2': df_inputgenes[1], 'name': 'Nonsense vs Silent'}, 
        {'df1': df_inputgenes[3], 'df2': df_inputgenes[1], 'name': 'No mutation vs Silent'}, 
                   ]
    for comp in comparisons: 
        if not comp['df1'].empty and not comp['df2'].empty: 
            U1, p = mannwhitneyu(comp['df1']['LFC'], comp['df2']['LFC'], method="asymptotic")
            comp['U1'] = U1
            comp['p'] = p
            print(f"{comp['name']}: {U1} {p}")

    df_InputGene_edits_list = pd.concat(df_inputgenes)
    df_InputGene_edits_list = df_InputGene_edits_list.reset_index(drop=True)

    return df_InputGene_edits_list, comparisons

def violin_plot(
        df_InputGene_edits_list, 
        edits_filedir, screen_name, input_gene, 
        directional=False, 
): 
    """
    Description
        Graph a violin plot of LFC distribution by category

    Params
        edits_filedir: Path, required
            Path to working directory
        screen_name: str, required
            the name of the input screen parsed form the input file
        input_gene: str, required
            the name of the input gene
        directional: bool, optional, default is False
            Whether or not to include bidirectional data

    Returns
        means: list of floats
            List of means for each mutation category
        stds: list of floats
            List of standard deviations for each mutation category
        medians: list of floats
            List of medians for each mutation category
    """

    fig, ax = plt.subplots()

    means = df_InputGene_edits_list.groupby('muttype')['LFC'].mean()
    stds = df_InputGene_edits_list.groupby('muttype')['LFC'].std()
    medians = df_InputGene_edits_list.groupby('muttype')['LFC'].median()

    if directional: 
        sns.violinplot(data=df_InputGene_edits_list, x="LFC", y="muttype", 
                       hue="LFC_direction", inner=None).set(title=screen_name)
        plotname = edits_filedir / f"plots/{input_gene}_{screen_name}_LFC_dist_muttype_bidirectional.pdf"
    else: 
        sns.violinplot(data=df_InputGene_edits_list, x="LFC", y="muttype", 
                       inner=None).set(title=screen_name)
        plotname = edits_filedir / f"plots/{input_gene}_{screen_name}_LFC_dist_muttype.pdf"
        
    plt.axvline(df_InputGene_edits_list["LFC"].mean(), c="gray", linestyle="dashed")
    plt.setp(ax.collections, alpha=.4)
    plt.scatter(y=range(len(means)), x=means, c="violet", alpha=.9)

    plt.savefig(plotname, dpi=300)

    return means, stds, medians
