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
comparisons = [
    ('Nonsense', 'Missense'), 
    ('Nonsense', 'No_Mutation'), 
    ('Nonsense', 'Silent'), 
    ('No_Mutation', 'Silent'), 
]

def mann_whitney_test(
    edits_filedir, screen_names, input_gene, 
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
    muts_dicts_list = [{} for _ in mut_categories_unspaced]
    for mut, mut_dict in zip(mut_categories_unspaced, muts_dicts_list): 
        list_mut = []
        for screen_name in screen_names: 
            edits_filename = f"screendata/{input_gene}_{screen_name}_{mut}_edits_list.tsv"
            df_temp = pd.read_csv(edits_filedir / edits_filename, sep = '\t')
            list_mut.extend(df_temp['LFC'].tolist())
            del df_temp
        mut_dict['LFC'] = list_mut
        mut_dict['muttype'] = mut
    muts_dicts = dict(map(lambda i, j : (i, j) , mut_categories_unspaced, muts_dicts_list))

    # MANN WHITNEY TEST #
    mannwhiteney_results = {}
    for comp1, comp2 in comparisons: 
        if len(muts_dicts[comp1]['LFC']) > 0 and len(muts_dicts[comp2]['LFC']) > 0: 
            U1, p = mannwhitneyu(muts_dicts[comp1]['LFC'], muts_dicts[comp2]['LFC'], method="asymptotic")
            mannwhiteney_results[f'{comp1} vs {comp2}'] = {'U1': U1, 'p': p}

    df_muts = pd.concat([pd.DataFrame(d) for d in muts_dicts_list]).reset_index(drop=True)
    return df_muts, mannwhiteney_results

def violin_plot(
        df_muts, edits_filedir, input_gene, screen_name, 
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
    fig, ax = plt.subplots(1, 2, figsize=(12,6))

    means = df_muts.groupby('muttype')['LFC'].mean()
    stds = df_muts.groupby('muttype')['LFC'].std()
    medians = df_muts.groupby('muttype')['LFC'].median()

    sns.violinplot(ax=ax[0], data=df_muts, x="LFC", y="muttype", 
                    inner=None).set(title='LFC by Mutation Type')
    sns.violinplot(ax=ax[1], data=df_muts, x="LFC", y="muttype", 
                    hue="LFC_direction", inner=None).set(title='LFC by Mutation Type')
    plt.axvline(df_muts["LFC"].mean(), c="gray", linestyle="dashed")
    plt.scatter(y=range(len(means)), x=means, c="violet", alpha=.9)
    plt.suptitle(screen_name)

    plt.tight_layout()
    plotname = edits_filedir / f"plots/{input_gene}_LFC_dist_muttype.pdf"
    plt.savefig(plotname, dpi=300)

    return means, stds, medians

def counts_by_gene(
    df_inputs, 
    gene_col, mut_col, 
    edits_filedir, title, 
): 
    # FIND HOW MANY PLOTS NEEDED #
    unique_genes = []
    for df_input in df_inputs: 
        unique = df_input[gene_col].unique().tolist()
        unique_genes = list(set(unique_genes+unique))
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))

    # MUTATION COUNTS ACROSS SCREENS BY GENE #
    df_mutation_counts = pd.DataFrame(columns = ['Gene'] + mut_categories_spaced, index=[0])
    for idx, gene in enumerate(unique_genes): 
        res = [gene] + [0 for _ in range(len(mut_categories_spaced))]
        for df_input in df_inputs:  
            df_current_gene = df_input.loc[df_input[gene_col] == gene,]
            for j, mutcat in enumerate(mut_categories_spaced): 
                res[j+1] += len(df_current_gene.loc[df_current_gene[mut_col] == mutcat, ])
        df_mutation_counts.loc[idx] = res

    df_plot = df_mutation_counts.melt("Gene", var_name="Mut Type", value_name="Count")
    # BARPLOT #
    sns.set_style("darkgrid")
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))
    ax = sns.catplot(data=df_plot, kind="bar", col_wrap=plot_dim, 
                     x="Mut Type", y="Count", hue="Mut Type", col="Gene", sharex = False)
    for ax in ax.axes.flat:
        for idx in range(len(ax.containers)):
            ax.bar_label(ax.containers[idx])
    plt.suptitle(title)

    # SAVE BARPLOT #
    plotname = f"plots/barplot_count_by_muttype.pdf"
    plt.savefig(edits_filedir / plotname, dpi=300)
    del df_mutation_counts, df_plot

def violin_by_gene(
    df_inputs, 
    gene_col, mut_col, val_col, 
    edits_filedir, title
): 
    # FIND HOW MANY PLOTS NEEDED #
    unique_genes = []
    for df_input in df_inputs: 
        unique = df_input[gene_col].unique().tolist()
        unique_genes = list(set(unique_genes+unique))
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))
    
    # VIOLIN PLOT SETUP #
    fig, axes = plt.subplots(nrows=plot_dim, ncols=plot_dim, sharex=False, sharey=True, 
                             figsize=(19,17), gridspec_kw={'hspace':0.3, 'wspace':0.1})

    for idx, current_gene in enumerate(unique_genes): 
        df_gene = pd.DataFrame()
        # AGGREGATE OVER EVERY SCREEN #
        for df_input in df_inputs: 
            df_current_gene = df_input.loc[df_input[gene_col] == current_gene,]

            # ENSURE ALL MUTATION TYPES PRESENT #
            for mutcat in mut_categories_spaced: 
                if mutcat not in df_current_gene[mut_col].unique():
                    df_temp = pd.DataFrame({val_col: [np.nan], mut_col: [mutcat]})
                    df_current_gene = pd.concat([df_current_gene, df_temp], ignore_index=True)
            df_gene = pd.concat([df_gene, df_current_gene], ignore_index=True)
        del df_current_gene

        # CALC MEAN STD #
        Means = df_gene.groupby(mut_col)[val_col].mean()
        # STDs = df_gene.groupby(mut_col)[val_col].std()
        # print(f"{current_gene}: Mean + STD")
        # for mutcat in mut_categories_spaced: 
        #     print(f"{mutcat}: {round(Means[mutcat], 3)} + {round(STDs[mutcat], 3)}")

        # PLOT VIOLIN #
        ax = axes.flatten()[idx]
        df_gene.loc[:, mut_col] = pd.Categorical(df_gene[mut_col], categories=mut_categories_spaced)
        df_gene = df_gene.sort_values(by=[mut_col]).reset_index(drop=True)
        sns.violinplot(ax=ax, data=df_gene, x=val_col, y=mut_col, 
                       inner=None, hue=mut_col).set(title=current_gene) # VIOLIN PLOTS #
        ax.axvline(df_gene[val_col].mean(), c="gray", linestyle="dashed")
        ax.scatter(y=range(len(Means)), x=Means, c="violet", alpha=.9) # MEANS #
        del df_gene
    plt.suptitle(title)

    # SAVE VIOLIN #
    plotname = f"plots/violinplot_LFC_by_muttype.pdf"
    plt.savefig(edits_filedir / plotname, dpi=300)
