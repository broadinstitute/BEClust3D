"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.1

"""

import pandas as pd
import re
from pathlib import Path

def parse_edit_helper(mut_type, this_edit): 
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
    match mut_type: 
        case "Missense": 
            assert edit_refAA != edit_altAA
        case "Silent": 
            assert edit_refAA == edit_altAA
        case "Nonsense": 
            assert edit_altAA in ['Ter', 'STOP', '*'] ### some possible ways to decribe stop codon

    return edit_refAA, edit_altAA, edit_pos

def parse_base_editing_results(df_InputGene, workdir, 
                               input_gene, input_screen, 
                               mut_col='Mutation category', val_col='logFC', 
                               gene_col='Target Gene Symbol', edits_col='Amino Acid Edits',
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
            each dataframe with column headers edit, human_pos, refAA, altAA, LFC
    """

    ### figure out creating directories

    # df_InputGene_missense, df_InputGene_silent, df_InputGene_nonsense, df_InputGene_intron, df_InputGene_UTR
    mut_categories = ['Missense', 'Silent', 'Nonsense', 'Intron', 'UTR'] # based on the standard output from BEAGLE
    mut_dfs = [pd.DataFrame() for _ in mut_categories]

    edits_filedir = Path(workdir + input_gene)

    for mut, df in zip(mut_categories, mut_dfs): 
        # preprocess by taking the subset of the original df with the indicated mutation
        df = df_InputGene.loc[df_InputGene[mut_col] == mut, ] ### == or .contains
        df = df.reset_index(drop=True)
        print(f"Count of {mut} rows: " + str(len(df)))

        # only Missense, Silent, Nonsense have associated mutations to parse
        if mut not in ['Missense', 'Silent', 'Nonsense']: 
            continue

        # open a file to write in list of mutations and their LFC values
        edits_filename = edits_filedir / f"screendata/{input_gene}_{input_screen}_{mut}_edits_list.tsv"
        edits_file = open(edits_filename, "w")
        edits_file.write("edit\thuman_pos\trefAA\taltAA\tLFC\n")

        # iterate through each cell of each list of edits
        for i in range(0, len(df)):
            edits_all = df.at[i, edits_col].strip(',').strip(';')
            edit_val = round(df.at[i, val_col], 3)
            edits_list = edits_all.split(',')
            
            for j in range(0, len(edits_list)):
                # identify edit information and output them into .tsv
                this_edit = edits_list[j].strip()
                edit_refAA, edit_altAA, edit_pos = parse_edit_helper(mut, this_edit)
                
                line = this_edit + '\t' + str(edit_pos) + '\t' + edit_refAA + '\t' + edit_altAA + '\t' + str(edit_val) + '\n'
                edits_file.write(line)
                
        edits_file.close()

    # preprocess No mutations separately
    df_InputGene_nomutation = df_InputGene.loc[df_InputGene[mut_col].isna(), ]
    df_InputGene_nomutation = df_InputGene_nomutation.reset_index(drop=True)
    print("Count of No mutation rows: " + str(len(df_InputGene_nomutation)))
    # no mutation
    df_nomutation_new = pd.DataFrame()
    df_nomutation_new['gene'] = df_InputGene_nomutation[gene_col]
    df_nomutation_new['LFC'] = df_InputGene_nomutation[val_col]
    df_nomutation_new.to_csv(edits_filedir / f"screendata/{input_gene}_{input_screen}_No_mutation_edits_list.tsv", index=False);

    ### did not include splice site mutations

    mut_dfs = mut_dfs.append(df_InputGene_nomutation)
    return mut_dfs
