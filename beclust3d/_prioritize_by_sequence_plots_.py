"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Helper Functions for Prioritize by Sequence

"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def counts_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, edit_type, 
): 
    # PREP DATA #
    counts = df_struc_consvr[f'all_{edit_type}_edits'].str.count(';').fillna(0).astype(int)+1
    counts[df_struc_consvr[f'all_{edit_type}_edits'] == '-'] = 0

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax = sns.barplot(x=df_struc_consvr['unipos'], y=counts, 
                     color='steelblue', edgecolor='steelblue')
    ax.set_ylabel(f"Count of {edit_type} Mutations")
    ax.set_xlabel(f"unipos")
    ax.set_title(f"{input_gene} Count of {edit_type} Mutations {screen_name}")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    counts_filename = f"plots/{input_gene}_{screen_name}_num_{edit_type}_per_residue.pdf"
    plt.savefig(edits_filedir / counts_filename, dpi=300)

def scatterplot_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, 
    edit_type, function_type, input='', 
): 
    # PREP DATA #
    x_list = df_struc_consvr['unipos'].tolist()
    y_list = df_struc_consvr[f'{function_type}_{edit_type}_LFC{input}'].tolist()
    x_vals = [x for x, y in zip(x_list, y_list) if y!='-']
    y_vals = [float(y) for y in y_list if y!='-']

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.axhline(-1.0, c="red", linestyle="--")
    ax.axhline(1.0, c="blue", linestyle="--")
    ax.axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=ax, x=x_vals, y=y_vals, color='steelblue', edgecolor='steelblue')
    ax.set_ylabel(f"{edit_type} LFC{input} Score")
    ax.set_xlabel(f"unipos")
    ax.set_title(f'{input_gene} {edit_type} LFC{input} Score By Residue {screen_name}')
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    scatter_filename = f"plots/{input_gene}_{screen_name}_{edit_type}_lfc{input}_score_by_residue.pdf"
    plt.savefig(edits_filedir / scatter_filename, dpi=300)

def dual_scatterplot_by_residue(
    df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type='Missense', 
): 
    df_struc_consvr = df_struc_consvr[df_struc_consvr[f'mean_{edit_type}_LFC'] != '-']
    df_struc_consvr_pos = df_struc_consvr[df_struc_consvr[f'mean_{edit_type}_LFC'] > 0.0]
    df_struc_consvr_neg = df_struc_consvr[df_struc_consvr[f'mean_{edit_type}_LFC'] < 0.0]

    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 6))
    for ax in axs: 
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

    axs[0].axhline(-1.0, c="red", linestyle="--")
    axs[0].axhline(1.0, c="blue", linestyle="--")
    axs[0].axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=axs[0], data=df_struc_consvr_pos, x="unipos", 
                    y=f'mean_{edit_type}_LFC_Z', hue=f'mean_{edit_type}_LFC_plab', palette='tab10')
    axs[0].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    axs[0].set_title(f'Positive LFC Values')

    axs[1].axhline(-1.0, c="red", linestyle="--")
    axs[1].axhline(1.0, c="blue", linestyle="--")
    axs[1].axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=axs[1], data=df_struc_consvr_neg, x="unipos", 
                    y=f'mean_{edit_type}_LFC_Z', hue=f'mean_{edit_type}_LFC_plab', palette='tab10')
    axs[1].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    axs[1].set_title(f'Negative LFC Values')

    plt.subplots_adjust(wspace=0.1)
    plt.suptitle(f'{input_gene} LFC_Z Score {screen_name}')

    scatter_filename = f"plots/{input_gene}_{screen_name}_{edit_type}_lfc_z_score_by_residue_posneg.pdf"
    plt.savefig(edits_filedir / scatter_filename, dpi=300)

def dual_histogram_by_residue(
    df_struc_consvr, edits_filedir, input_gene, screen_name, edit_type='Missense', 
):  
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 6))
    for ax in axs: 
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

    df_struc_consvr = df_struc_consvr[df_struc_consvr[f'mean_{edit_type}_LFC'] != '-']
    df_struc_consvr_pos = df_struc_consvr[df_struc_consvr[f'mean_{edit_type}_LFC'] > 0.0]
    df_struc_consvr_neg = df_struc_consvr[df_struc_consvr[f'mean_{edit_type}_LFC'] < 0.0]
    plot1 = sns.histplot(ax=axs[0], data=df_struc_consvr_pos, x=f'mean_{edit_type}_LFC', hue=f'mean_{edit_type}_LFC_plab', bins=80, palette='tab10')
    plot2 = sns.histplot(ax=axs[1], data=df_struc_consvr_neg, x=f'mean_{edit_type}_LFC', hue=f'mean_{edit_type}_LFC_plab', bins=80, palette='tab10')

    # give corresponding titles
    plot1.set_title(f'Positive LFC Counts')
    plot2.set_title(f'Negative LFC Counts')
    plt.suptitle(f'{input_gene} Mean Missense LFC Counts {screen_name}')
    plt.subplots_adjust(wspace=0.1)

    hist_filename = f"plots/{input_gene}_{screen_name}_{edit_type}_lfc_score_by_bin_posneg.pdf"
    plt.savefig(edits_filedir / hist_filename, dpi=300)