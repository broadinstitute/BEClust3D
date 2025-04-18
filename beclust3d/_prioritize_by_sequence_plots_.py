"""
File: .py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-25
Description: Helper Functions for Prioritize by Sequence

"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

def counts_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, mut, 
): 
    # PREP DATA #
    counts = df_struc_consvr[f'all_{mut}_edits'].str.count(';').fillna(0).astype(int)+1
    counts[df_struc_consvr[f'all_{mut}_edits'] == '-'] = 0

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax = sns.barplot(x=df_struc_consvr['unipos'], y=counts, 
                     color='steelblue', edgecolor='steelblue')
    ax.set_ylabel(f"Count of {mut} Mutations")
    ax.set_xlabel(f"unipos")
    ax.set_title(f"{input_gene} Count of {mut} Mutations {screen_name}")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    counts_filename = f"plots/{input_gene}_{screen_name}_num_{mut}_per_res.pdf"
    plt.savefig(edits_filedir / counts_filename, dpi=300)
    # plt.close(fig)

def stdev_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, function_name, mut, yaxis=True,
): 
    # PREP DATA #
    xvals = df_struc_consvr['unipos']
    yvals = pd.to_numeric(df_struc_consvr[f'{function_name}_{mut}_LFC'], errors='coerce').fillna(0)
    stdevs = pd.to_numeric(df_struc_consvr[f'{function_name}_{mut}_LFC_stdev'], errors='coerce').fillna(0)
    xvals_filtered = xvals[stdevs != 0]
    if yaxis: yvals_filtered = yvals[stdevs != 0]
    else: yvals_filtered = [0]*len(xvals_filtered)
    stdevs_filtered = stdevs[stdevs != 0]

    # PLOT #
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.set_facecolor('#EBEBEB')
    [ax.spines[side].set_visible(False) for side in ax.spines]
    ax.grid(which='major', color='white', linewidth=0.5)
    ax.set_axisbelow(True)

    ax.errorbar(x=xvals_filtered, y=yvals_filtered, yerr=stdevs_filtered, 
                color='steelblue', ls=' ', marker='o', capsize=3, capthick=1, ecolor='black', 
                # fmt='-', color='steelblue', ecolor='steelblue'
                )
    ax.set_ylabel(f"Standard Deviations of {mut} Mutations")
    ax.set_xlabel(f"unipos")
    ax.set_title(f"{input_gene} Standard Dev of {mut} Mutations {screen_name}")
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    stdev_filename = f"plots/{input_gene}_{screen_name}_stdev_{mut}_per_res.pdf"
    plt.savefig(edits_filedir / stdev_filename, dpi=300)
    # plt.close(fig)
    
def scatterplot_by_residue(
    df_struc_consvr, 
    edits_filedir, input_gene, screen_name, 
    mut, function_type, input='', 
): 
    # PREP DATA #
    x_list = df_struc_consvr['unipos'].tolist()
    y_list = df_struc_consvr[f'{function_type}_{mut}_LFC{input}'].tolist()
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
    ax.set_ylabel(f"{mut} LFC{input} Score")
    ax.set_xlabel(f"unipos")
    ax.set_title(f'{input_gene} {mut} LFC{input} Score By Residue {screen_name}')
    plt.xticks(np.arange(0, len(df_struc_consvr), 50), rotation = 90)

    scatter_filename = f"plots/{input_gene}_{screen_name}_{mut}_lfc{input}_score_by_res.pdf"
    plt.savefig(edits_filedir / scatter_filename, dpi=300)
    # plt.close(fig)

def dual_scatterplot_by_residue(
    df_struc_consvr, edits_filedir, input_gene, screen_name, function_type, mut='Missense', 
): 
    df_struc_consvr = df_struc_consvr[df_struc_consvr[f'{function_type}_{mut}_LFC'] != '-']
    df_struc_consvr[f'{function_type}_{mut}_LFC'] = df_struc_consvr[f'{function_type}_{mut}_LFC'].astype(float)
    df_struc_consvr_pos = df_struc_consvr[df_struc_consvr[f'{function_type}_{mut}_LFC'] > 0.0]
    df_struc_consvr_neg = df_struc_consvr[df_struc_consvr[f'{function_type}_{mut}_LFC'] < 0.0]

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
                    y=f'{function_type}_{mut}_LFC_Z', hue=f'{function_type}_{mut}_LFC_plab', palette='tab10')
    axs[0].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    axs[0].set_title(f'Positive LFC Values')

    axs[1].axhline(-1.0, c="red", linestyle="--")
    axs[1].axhline(1.0, c="blue", linestyle="--")
    axs[1].axhline(0.0, c="gray", linestyle="--")
    sns.scatterplot(ax=axs[1], data=df_struc_consvr_neg, x="unipos", 
                    y=f'{function_type}_{mut}_LFC_Z', hue=f'{function_type}_{mut}_LFC_plab', palette='tab10')
    axs[1].legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0)
    axs[1].set_title(f'Negative LFC Values')

    plt.subplots_adjust(wspace=0.3)
    plt.suptitle(f'{input_gene} LFC_Z Score {screen_name}')

    scatter_filename = f"plots/{input_gene}_{screen_name}_{mut}_lfcz_scatter_by_bin_posneg.pdf"
    plt.savefig(edits_filedir / scatter_filename, dpi=300)
    # plt.close(fig)

def dual_histogram_by_residue(
    df_struc_consvr, edits_filedir, input_gene, screen_name, function_type, mut='Missense', 
):  
    fig, axs = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(12, 6))
    for ax in axs: 
        ax.set_facecolor('#EBEBEB')
        [ax.spines[side].set_visible(False) for side in ax.spines]
        ax.grid(which='major', color='white', linewidth=0.5)
        ax.set_axisbelow(True)

    df_struc_consvr = df_struc_consvr[df_struc_consvr[f'{function_type}_{mut}_LFC'] != '-']
    df_struc_consvr[f'{function_type}_{mut}_LFC'] = df_struc_consvr[f'{function_type}_{mut}_LFC'].astype(float)
    df_struc_consvr_pos = df_struc_consvr[df_struc_consvr[f'{function_type}_{mut}_LFC'] > 0.0]
    df_struc_consvr_neg = df_struc_consvr[df_struc_consvr[f'{function_type}_{mut}_LFC'] < 0.0]
    plot1 = sns.histplot(ax=axs[0], data=df_struc_consvr_pos, x=f'{function_type}_{mut}_LFC', hue=f'{function_type}_{mut}_LFC_plab', bins=80, palette='tab10')
    plot2 = sns.histplot(ax=axs[1], data=df_struc_consvr_neg, x=f'{function_type}_{mut}_LFC', hue=f'{function_type}_{mut}_LFC_plab', bins=80, palette='tab10')

    # give corresponding titles
    plot1.set_title(f'Positive LFC Counts')
    plot2.set_title(f'Negative LFC Counts')
    plt.suptitle(f'{input_gene} Mean Missense LFC Counts {screen_name}')
    plt.subplots_adjust(wspace=0.1)

    hist_filename = f"plots/{input_gene}_{screen_name}_{mut}_lfc_hist_by_bin_posneg.pdf"
    plt.savefig(edits_filedir / hist_filename, dpi=300)
    # plt.close(fig)
