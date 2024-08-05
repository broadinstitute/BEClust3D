"""
File: preprocess_be_results.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-08-02
Description: Translated from Notebook 3.1

"""

import pandas as pd
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math

mut_categories_spaced = ["Nonsense", "Splice Site", "Missense", "No Mutation", "Silent"]
mut_categories_unspaced = [mc.replace(' ' , '_') for mc in mut_categories_spaced]
comparisons = [(mut_categories_unspaced[0], mut_categories_unspaced[2]), 
               (mut_categories_unspaced[0], mut_categories_unspaced[3]), 
               (mut_categories_unspaced[0], mut_categories_unspaced[4]), 
               (mut_categories_unspaced[3], mut_categories_unspaced[4]), 
               ]

def mann_whitney_test(
    edits_filedir, screen_name, input_gene, 
): 
    """
    Description
        A helper function to run the Mann Whitney test on
        'Missense', 'Silent', 'Nonsense', 'No Mutation'
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
    df_inputgenes = [pd.DataFrame() for _ in mut_categories_unspaced]
    for mut, df in zip(mut_categories_unspaced, df_inputgenes): 
        edits_filename = f"screendata/{input_gene}_{screen_name}_{mut}_edits_list.tsv"
        print(edits_filename)
        df_temp = pd.read_csv(edits_filedir / edits_filename, sep = '\t')
        df['LFC'] = df_temp['LFC']
        df['muttype'] = mut
    df_dict = dict(map(lambda i, j : (i, j) , mut_categories_unspaced, df_inputgenes))

    # MANN WHITNEY TEST #
    mannwhiteney_results = {}
    for comp1, comp2 in comparisons: 
        if not df_dict[comp1].empty and not df_dict[comp2].empty: 
            U1, p = mannwhitneyu(df_dict[comp1]['LFC'], df_dict[comp2]['LFC'], method="asymptotic")
            mannwhiteney_results[f'{comp1} vs {comp2}'] = {'U1': U1, 'p': p}

    df_InputGene_edits_list = pd.concat(df_inputgenes).reset_index(drop=True)

    return df_InputGene_edits_list, mannwhiteney_results

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

def counts_violin_by_gene(
    df_rawinput, 
    gene_col, mut_col, val_col, 
    edits_filedir, screen_name, 
): 
    unique_genes = df_rawinput[gene_col].unique()

    # MUTATION COUNTS BY GENE #
    df_mutation_counts = pd.DataFrame(columns = ['Gene'] + mut_categories_spaced, index=[0])
    for idx, gene in enumerate(unique_genes):
        df_current_gene = df_rawinput.loc[df_rawinput[gene_col] == gene,]
        res = [gene]
        for mutcat in mut_categories_spaced: 
            res.append(len(df_current_gene.loc[df_current_gene[mut_col] == mutcat, ]))
        df_mutation_counts.loc[idx] = res

    # BARPLOT #
    sns.set_style("darkgrid")
    plt.tight_layout()
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))
    ax = sns.catplot(data=df_mutation_counts.melt("Gene", var_name="Mutation Type", value_name="Count"), 
                     x="Mutation Type", y="Count", kind="bar", col_wrap=plot_dim, hue="Mutation Type", 
                     col="Gene", sharex = False )
    for ax in ax.axes.flat:
        for idx in range(len(ax.containers)):
            ax.bar_label(ax.containers[idx])
    del df_mutation_counts

    # SAVE BARPLOT #
    plotname = f"plots/{screen_name}_muttype_count.pdf"
    plt.savefig(edits_filedir / plotname, dpi=500)

    plot_dim = math.ceil(math.sqrt(len(unique_genes)))
    fig, axes = plt.subplots(nrows=plot_dim, ncols=plot_dim, sharex=False, sharey=True, 
                             figsize=(19,17), gridspec_kw={'hspace': 0.3, 'wspace': 0.1})

    for idx, current_gene in enumerate(unique_genes): 
        df_current_gene = df_rawinput.loc[df_rawinput[gene_col] == gene,]

        # ENSURE ALL MUTATION TYPES PRESENT #
        for mutcat in mut_categories_spaced: 
            if mutcat not in df_current_gene[mut_col].unique():
                df_current_gene = pd.concat([df_current_gene, 
                                             pd.DataFrame({val_col: [np.nan], mut_col: [mutcat]})], ignore_index=True)

        Means = df_current_gene.groupby(mut_col)[val_col].mean()
        STDs = df_current_gene.groupby(mut_col)[val_col].std()
        Medians = df_current_gene.groupby(mut_col)[val_col].median()

        print(f"{current_gene}: Mean + STD")
        for mutcat in mut_categories_spaced: 
            print(f"{mutcat}: {round(Means[mutcat], 3)} + {round(STDs[mutcat], 3)}") # Access by label for robustness

        # PLOT VIOLIN #
        ax = axes.flatten()[idx]
        df_current_gene.loc[:, mut_col] = pd.Categorical(df_current_gene[mut_col], categories=mut_categories_spaced)
        df_current_gene = df_current_gene.sort_values(by=[mut_col]).reset_index(drop=True)
        sns.violinplot(ax=ax, data=df_current_gene, x=val_col, y=mut_col, 
                       inner=None, hue=mut_col).set(title=current_gene)
        ax.axvline(df_current_gene[val_col].mean(), c="gray", linestyle="dashed")
        ax.scatter(y=range(len(Means)), x=Means, c="violet", alpha=.9)

    # SAVE VIOLIN #
    plotname = f"plots/{screen_name}_muttype_LFC_dist.pdf"
    plt.savefig(edits_filedir / plotname, dpi=500)
