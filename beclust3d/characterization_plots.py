"""
File: characterization_plots.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-04
Description: Plotting functions for hit characterization

"""

import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.lines import Line2D

def lfc_vs_lfc3d_scatterplot(
        df_lfc3d_dis, df_nonaggr_lfc3d, workdir, input_gene, screen_name, plot_name, lfc3d_hit_threshold=0.05
):
    """
    Description
        Generate LFC vs LFC3D scatter plot

    Params
        df_lfc3d_dis: DataFrame
            LFC3D distribution DataFrame
        df_nonaggr_lfc3d: DataFrame
            Non-aggregated LFC3D DataFrame
        workdir: str
            Working directory
        input_gene: str
            Input gene name
        screen_name: str
            Screen name
        plot_name: str
            Plot name to be saved in plots folder
        lfc3d_hit_threshold: float
            LFC3D hit threshold for determining significance
    """

    edits_filedir = Path(workdir)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    if not os.path.exists(edits_filedir / 'plots'):
        os.mkdir(edits_filedir / 'plots')


    # Load LFC and LFC3D scores and distributions
    df_lfc3d_dis.rename(columns={
        f"{screen_name}_LFC": "LFC",
        f"{screen_name}_LFC3D": "LFC3D",
        f"{screen_name}_LFC3D_dis": "LFC3D_dis"
    }, inplace=True)

    # Get LFC3D significance labels
    df_nonaggr_lfc3d.rename(columns={
        f'{screen_name}_LFC3D_neg_psig': "LFC3D_neg_psig",
        f'{screen_name}_LFC3D_pos_psig': "LFC3D_pos_psig"
    }, inplace=True)
    df_lfc3d_psig = df_nonaggr_lfc3d[['unipos', 'LFC3D_neg_psig', 'LFC3D_pos_psig']]
    df_lfc_lfc3d = pd.merge(df_lfc3d_dis, df_lfc3d_psig, on='unipos', how='left')

    # Assign p-significance label for hue coloring
    psig_dict = {'above': f'p>={lfc3d_hit_threshold}', 'below': f'p<{lfc3d_hit_threshold}'}

    def assign_psig_label(row):
        if row['LFC3D_neg_psig'] == psig_dict['above'] and row['LFC3D_pos_psig'] == psig_dict['above']:
            return 'not hit'
        elif row['LFC3D_neg_psig'] == psig_dict['above'] and row['LFC3D_pos_psig'] == psig_dict['below']:
            return 'positive hit'
        elif row['LFC3D_neg_psig'] == psig_dict['below'] and row['LFC3D_pos_psig'] == psig_dict['above']:
            return 'negative hit'
        elif row['LFC3D_neg_psig'] == psig_dict['below'] and row['LFC3D_pos_psig'] == psig_dict['below']:
            return 'pos/neg hit'
        return None
    df_lfc_lfc3d['psig_label'] = df_lfc_lfc3d.apply(assign_psig_label, axis=1)

    # Remove dashes from table and replace with 0
    df_lfc_lfc3d['LFC'] = df_lfc_lfc3d['LFC'].replace('-', 0.0).astype(float)
    df_lfc_lfc3d['LFC3D'] = df_lfc_lfc3d['LFC3D'].replace('-', 0.0).astype(float)

    y_min = df_lfc_lfc3d['LFC3D'].min()
    x_min = df_lfc_lfc3d['LFC'].min()
    df_lfc_lfc3d['LFC'] = df_lfc_lfc3d['LFC'].replace(0.0, x_min-1).astype(float)
    # print(input_gene)
    # display(df_combined.head(7))

    # Hit Type Colors
    custom_palette = {
        'not hit': 'grey',
        'positive hit': 'blue',
        'negative hit': 'red',
        'pos/neg hit': 'magenta'
    }
    # Scatter plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=df_lfc_lfc3d, x='LFC', y='LFC3D', hue="psig_label", palette=custom_palette)
    plt.axhline(y_min, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(x_min, color="gray", linestyle="--", linewidth=0.8)
    plt.title(f"{input_gene} LFC vs LFC3D Scatter Plot")
    plt.xlabel(f"{screen_name} (LFC)")
    plt.ylabel(f"{screen_name} (LFC3D)")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.savefig(plot_name, dpi=300)

def RSA_vs_pLDDT_barplot(df_filtered, gene_name, plot_name):
    """
    Description
        Generate RSA vs pLDDT barplot
    
    Params
        df_filtered: DataFrame
            Filtered DataFrame with relevant colummns (bfactor_pLDDT, RSA, LFC3D_wght, dir)
        gene_name: str
            Gene name
        plot_name: str
            Plot name to be saved in plots folder
    """

    color_map = {'NEG': 'darkred', 'POS': 'darkblue'}
    colors = df_filtered['dir'].map(color_map)

    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(
        df_filtered['bfactor_pLDDT'],
        df_filtered['RSA'],
        s=df_filtered['LFC3D_wght'],
        c=colors, alpha=0.7)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='POS',
            markerfacecolor='darkred', markersize=10),
        Line2D([0], [0], marker='o', color='w', label='NEG',
            markerfacecolor='darkblue', markersize=10)
    ]

    sizes = [5, 50, 95]
    for size in sizes:
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', label=f'Size {size}',
                   markerfacecolor='gray', markersize=np.sqrt(size)) )

    plt.legend(handles=legend_elements, title="Legend")
    plt.xlabel('pLDDT')
    plt.ylabel('RSA')
    plt.title(f"{gene_name} RSA vs. pLDDT Scatterplot")
    plt.savefig(plot_name, dpi=300)

def hits_vs_feature_barplot(df_filtered, xcolumn, xname, gene_name, plot_name, type='count'):
    """
    Description
        Generate hit count barplot for specified feature (ex. RSA, pLDDT, Domain, etc.)
    
    Params
        df_filtered: DataFrame
            Filtered DataFrame with relevant colummns (domain, pLDDT_dis, exposure, SS3, etc...)
        xcolumn: str
            Column name for the feature
        xname: str
            Feature name
        gene_name: str
            Gene name
        plot_name: str
            Plot name to be saved in plots folder
        type: str, optional
            Type of plot: 'count' for direct hit counts, 'fraction' for fraction percentages (default: 'count')
    """
    plt.figure(figsize=(10, 6))

    if type == 'fraction':
        # Calculate total counts of POS and NEG across the entire dataset
        total_pos = df_filtered[df_filtered['dir'] == 'POS'].shape[0]
        total_neg = df_filtered[df_filtered['dir'] == 'NEG'].shape[0]

        # Create a DataFrame with the original counts for each xcolumn category and dir
        count_data = df_filtered.groupby([xcolumn, 'dir']).size().unstack(fill_value=0).reset_index()

        # Calculate fractions for POS and NEG independently (normalize so each sums to 100%)
        if total_pos > 0:
            count_data['POS'] = (count_data['POS'] / total_pos) * 100
        if total_neg > 0:
            count_data['NEG'] = (count_data['NEG'] / total_neg) * 100

        # Melt the DataFrame to a long format for plotting
        fraction_data_melted = count_data.melt(id_vars=xcolumn, value_vars=['POS', 'NEG'], 
                                               var_name='dir', value_name='fraction')

        # Plot the fractions
        sns.barplot(data=fraction_data_melted, x=xcolumn, y='fraction', hue='dir', 
                    palette={'NEG': 'darkred', 'POS': 'darkblue'})
        plt.ylabel('Fraction of Hits (%)')
    else:
        # Count plot
        sns.countplot(data=df_filtered, x=xcolumn, hue='dir', 
                      palette={'NEG': 'darkred', 'POS': 'darkblue'})
        plt.ylabel('Count of Hits')

    plt.xlabel(xname)
    plt.title(f"{gene_name} {xname} Hit Count Barplot")
    plt.legend(title='dir')
    plt.savefig(plot_name)