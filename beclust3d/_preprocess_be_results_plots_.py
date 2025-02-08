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
        Run the Mann Whitney test on 'Missense', 'Silent', 'Nonsense', 'No Mutation'
    """

    muts_dicts_list = [{} for _ in mut_categories_unspaced]
    for mut, mut_dict in zip(mut_categories_unspaced, muts_dicts_list): 
        list_mut = []
        for screen_name in screen_names: 
            edits_filename = f"screendata/{input_gene}_{screen_name}_{mut}.tsv"
            try:
                df_temp = pd.read_csv(edits_filedir / edits_filename, sep='\t')
                list_mut.extend(df_temp['LFC'].tolist())
                del df_temp
            except FileNotFoundError:
                list_mut.extend([])
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
    """

    plt.rcParams.update({'font.size': 10})
    fig, ax = plt.subplots(1, 2, figsize=(12,6))

    means = df_muts.groupby('muttype')['LFC'].mean()
    stds = df_muts.groupby('muttype')['LFC'].std()
    medians = df_muts.groupby('muttype')['LFC'].median()
    sorted_muttypes = sorted(means.sort_values().index.tolist())

    sns.violinplot(ax=ax[0], data=df_muts, x="LFC", y="muttype", 
                   order=sorted_muttypes, inner=None).set(title='LFC by Mutation Type')
    sns.violinplot(ax=ax[1], data=df_muts, x="LFC", y="muttype", 
                   order=sorted_muttypes, hue="LFC_direction", inner=None).set(title='LFC by Mutation Type')
    plt.axvline(df_muts["LFC"].mean(), c="gray", linestyle="dashed")
    plt.scatter(y=range(len(means)), x=means, c="violet", alpha=.9)
    plt.suptitle(f'{screen_name}_{input_gene}')

    plt.tight_layout()
    plotname = edits_filedir / f"plots/{screen_name}_{input_gene}_LFC_dist_by_muttype.pdf"
    plt.savefig(plotname, dpi=300)

    return means, stds, medians


def counts_by_gene(
    df_inputs, 
    gene_col, mut_col, 
    edits_filedir, title, 
): 
    """
    Description
        Graph a bar plot of counts by category per gene
    """

    mut_categories_spaced_sort = sorted(mut_categories_spaced)

    # Compute unique genes efficiently
    all_genes = pd.concat([df[gene_col] for df in df_inputs], ignore_index=True)
    unique_genes = sorted(all_genes.unique().tolist())
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))

    # Compute mutation counts for each gene and mutation category
    df_all = pd.concat(df_inputs, ignore_index=True)
    df_mutation_counts = (
        df_all.groupby([gene_col, mut_col])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=mut_categories_spaced_sort, fill_value=0)
        .reset_index()
    )
    df_mutation_counts.columns = ['Gene'] + mut_categories_spaced_sort

    # # MUTATION COUNTS ACROSS SCREENS BY GENE #
    # df_mutation_counts = pd.DataFrame(columns = ['Gene'] + mut_categories_spaced_sort, index=[0])
    # for idx, gene in enumerate(unique_genes): 
    #     res = [gene] + [0 for _ in range(len(mut_categories_spaced_sort))]
    #     for df_input in df_inputs:  
    #         df_current_gene = df_input.loc[df_input[gene_col] == gene,]
    #         for j, mutcat in enumerate(mut_categories_spaced_sort): 
    #             res[j+1] += len(df_current_gene.loc[df_current_gene[mut_col] == mutcat, ])
    #     df_mutation_counts.loc[idx] = res

    # BARPLOT #
    df_plot = df_mutation_counts.melt("Gene", var_name="Mut Type", value_name="Count")
    sns.set_style("darkgrid")
    ax = sns.catplot(data=df_plot, kind="bar", col_wrap=plot_dim, 
                    x="Mut Type", y="Count", hue="Mut Type", col="Gene", sharex=False)

    for ax in ax.axes.flat:
        for container in ax.containers:
            ax.bar_label(container)

    plt.subplots_adjust(top=0.9)
    plt.suptitle(title)

    # SAVE BARPLOT #
    plotname = f"plots/{title}_barplot_count_by_muttype.pdf"
    plt.savefig(edits_filedir / plotname, dpi=300)
    del df_mutation_counts, df_plot

def violin_by_gene(
    df_inputs, 
    gene_col, mut_col, val_col, 
    edits_filedir, title
): 
    """
    Description
        Graph a violin plot of LFC values by category per gene
    """
    
    mut_categories_spaced_sort = sorted(mut_categories_spaced)
    # FIND HOW MANY PLOTS NEEDED #
    unique_genes = set()
    for df_input in df_inputs: 
        unique_genes.update(df_input[gene_col].unique().tolist())
    unique_genes = sorted(unique_genes)
    plot_dim = math.ceil(math.sqrt(len(unique_genes)))

    # VIOLIN PLOT SETUP #
    plt.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(nrows=plot_dim, ncols=plot_dim, sharex=False, sharey=True, 
                             figsize=(19,17), gridspec_kw={'hspace':0.3, 'wspace':0.1})
    axes = axes.flatten()

    legend_handles, legend_labels = None, None

    for idx, current_gene in enumerate(unique_genes): 
        df_gene = pd.DataFrame()
        # AGGREGATE OVER EVERY SCREEN #
        for df_input in df_inputs: 
            df_current_gene = df_input.loc[df_input[gene_col] == current_gene,]

            # ENSURE ALL MUTATION TYPES PRESENT #
            for mutcat in mut_categories_spaced_sort: 
                if mutcat not in df_current_gene[mut_col].unique():
                    df_temp = pd.DataFrame({val_col: [np.nan], mut_col: [mutcat]})
                    df_current_gene = pd.concat([df_current_gene, df_temp], ignore_index=True)
            df_gene = pd.concat([df_gene, df_current_gene], ignore_index=True)
        del df_current_gene

        # CALC MEAN STD #
        filtered_df_gene = df_gene[df_gene[mut_col].isin(mut_categories_spaced_sort)]
        Means = filtered_df_gene.groupby(mut_col)[val_col].mean()

        # PLOT VIOLIN #
        df_gene.loc[:, mut_col] = pd.Categorical(df_gene[mut_col], categories=mut_categories_spaced_sort)
        df_gene = df_gene.sort_values(by=[mut_col]).reset_index(drop=True)
        violin = sns.violinplot(ax=axes[idx], data=df_gene, x=val_col, y=mut_col, 
                                inner=None, hue=mut_col)
        axes[idx].set_title(current_gene)

        # Extract legend from the first plot (only once)
        if legend_handles is None and legend_labels is None:
            legend_handles, legend_labels = axes[idx].get_legend_handles_labels()

        # Remove legend from the individual plot
        axes[idx].get_legend().remove()
        axes[idx].axvline(df_gene[val_col].mean(), c="gray", linestyle="dashed")
        axes[idx].scatter(y=range(len(Means)), x=Means, c="violet", alpha=.9) # MEANS #
        del df_gene
    
    # REMOVE UNUSSED AXES
    for j in range(idx + 1, len(axes)):
        axes[j].set_visible(False)
    
    plt.suptitle(title)
    plt.subplots_adjust(top=0.9, wspace=0.1)

    # SAVE VIOLIN #
    plotname = f"plots/{title}_violinplot_LFC_by_muttype.pdf"
    plt.savefig(edits_filedir / plotname, dpi=300)

    # CREATE SEPARATE LEGEND #
    if legend_handles and legend_labels:
        legend_fig = plt.figure(figsize=(4, 2))  # Adjust size as needed
        legend_ax = legend_fig.add_subplot(111)
        legend_ax.axis('off')  # Hide axes
        legend_ax.legend(legend_handles, legend_labels, loc='center')
        
        legend_path = f"plots/{title}_legend.pdf"
        legend_fig.savefig(edits_filedir / legend_path, dpi=300)
