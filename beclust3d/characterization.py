"""
File: characterization.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-23
Description: Characterization of hits and enrichment tests

"""

import pandas as pd
from pathlib import Path
import os
import subprocess
import shutil
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.stats import fisher_exact

def add_annotation(df_coord_struc, workdir, input_gene, input_uniprot, structureid):
    # Download UniProt data
    ffile = f"{input_uniprot}.txt"
    target_path = os.path.join("drive/MyDrive/BEPipeline", input_gene)
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    subprocess.run(["wget", f"https://rest.uniprot.org/uniprotkb/{ffile}"])

    downloaded_file = f"./{ffile}"  # File in Colab's local environment
    if os.path.exists(downloaded_file):
        target_file = os.path.join(target_path, ffile)
        if not os.path.exists(target_file): # check if the destination file already exists
            shutil.move(downloaded_file, target_path)
            print(f"{downloaded_file} successfully moved to {target_path}")
        else:
            print(f"{target_file} already exists. Skipping move.")
    else:
        raise FileNotFoundError(f"File {downloaded_file} not found after download!")

    df_META = pd.DataFrame()
    df_META["unipos"] = df_coord_struc["unipos"]
    df_META["unires"] = df_coord_struc["unires"]

    uniprot_source_filename = f"{workdir}/{input_uniprot}.txt"
    uniprot_source_file = open(uniprot_source_filename, "r")
    uniprot_source_lines = uniprot_source_file.readlines()
    uniprot_source_file.close()

    uniprot_processed_lines = [idx for idx in uniprot_source_lines if idx[0:2] == "FT"]
    uniprot_processed_lines_joined = "\n".join(uniprot_processed_lines)

    uniprot_processed_filename = f"{workdir}/{input_uniprot}_processed.uniprot"
    uniprot_processed_file = open(uniprot_processed_filename, "w")
    uniprot_processed_file.writelines(uniprot_processed_lines)
    uniprot_processed_file.close()

    arr_domain = []
    for i in range(0,len(df_META)):
        unipos = df_META.at[i, 'unipos']
        uniaa = df_META.at[i, 'unires']

        # Domains
        ## domain (code: DOMAIN)
        domain = '-'
        for j in range(0, len(uniprot_processed_lines)):
            featureline = uniprot_processed_lines[j].strip()
            if "DOMAIN" in featureline:
                positionline = featureline.strip()
                positioninfo = positionline[21:len(positionline)]
                positiontokens = positioninfo.split('..')

                notesline = uniprot_processed_lines[j+1].strip()
                notesinfo = notesline[21:len(notesline)]
                notestoken = notesinfo.split('=')
                notes = notestoken[1].strip()
                notes = notes[1:len(notes)-1]

                if len(positiontokens) == 1:
                    position = int(positiontokens[0])
                    if unipos == position:
                        domain = str(positioninfo) + ":" + str(notes)
                elif len(positiontokens) == 2:
                    start_position = int(positiontokens[0])
                    end_position = int(positiontokens[1])
                    if (unipos >= start_position) and (unipos <= end_position):
                        domain = str(positioninfo) + ":" + str(notes)
                elif len(positiontokens) == 0:
                    print("feature position can't be empty!")
                    break;
                elif len(positiontokens) > 2:
                    print("feature positions can't be more than two!")
                    break;
        arr_domain.append(domain)

    df_META['domain'] = arr_domain

    df_META.to_csv(f"{workdir}/{input_gene}_{input_uniprot}_Uniprotf.tsv", sep = '\t', index=False)
    df_merged = pd.merge(df_coord_struc, df_META, on=['unipos', 'unires'], how='left')
    # Resave coord_struc_features with added domains
    df_merged.to_csv(f"{workdir}/{structureid}_coord_struc_features.tsv", sep = '\t', index=False)
    return df_merged

# Test 1: High pLDDT vs. Low pLDDT
def high_vs_low_pLDDT(data, pthr, domain_labels):
    below = data[data['mean_Missense_LFC_Z'] < pthr]
    above = data[data['mean_Missense_LFC_Z'] >= pthr]
    table = np.array([
        [len(below[below['pLDDT_dis'].isin(domain_labels['in'])]),
            len(below[below['pLDDT_dis'].isin(domain_labels['out'])])],
        [len(above[above['pLDDT_dis'].isin(domain_labels['in'])]),
            len(above[above['pLDDT_dis'].isin(domain_labels['out'])])]
    ])
    ftest = fisher_exact(table, alternative='two-sided')
    odds_ratio = stats.contingency.odds_ratio(table)
    ci = odds_ratio.confidence_interval(confidence_level=0.95)
    return ftest, odds_ratio, ci

 # Test 2: In-Domain vs. Outside-Domain
def inside_vs_outside_domain(data, pthr, domain_labels):
    below = data[data['mean_Missense_LFC_Z'] < pthr].reset_index(drop=True)
    above = data[data['mean_Missense_LFC_Z'] >= pthr].reset_index(drop=True)
    below_in_domain = below[below['domain'] != '-']
    below_out_domain = below[below['domain'] == '-']
    above_in_domain = above[above['domain'] != '-']
    above_out_domain = above[above['domain'] == '-']
    table = np.array([
        [len(below_in_domain), len(below_out_domain)],
        [len(above_in_domain), len(above_out_domain)]
    ])
    ftest = fisher_exact(table, alternative='two-sided')
    odds_ratio = stats.contingency.odds_ratio(table)
    ci = odds_ratio.confidence_interval(confidence_level=0.95)
    return ftest, odds_ratio, ci

def plot_enrichment_tests(results, pthr, input_gene):
    fig, ax = plt.subplots(figsize=(10, 6))
    y_positions = [1, 2, 3, 4]  # Fixed y-axis positions for the four lines

    for i, result in enumerate(results):
        odds_ratio = result['odds_ratio']
        ci = result['ci']
        y = y_positions[i]
        # color = 'blue' if result['lfc_type'] == 'Positive LFC' else 'red'
        color = 'red' if i % 2 == 0 else 'blue'
        if np.isnan(odds_ratio) or np.isinf(odds_ratio):
            # Placeholder for NaN or infinite odds ratio
            x_min, x_max = ax.get_xlim()
            x_mid = (x_min + x_max) / 2
            ax.plot([x_mid], [y], 'o', color="grey", markerfacecolor='none', linestyle='None')
            continue

        # Determine styling
        is_significant = result['p_value'] <= pthr
        marker_style = 'o' if is_significant else 'o'
        marker_fill = color if is_significant else 'none'
        line_style = '-' if is_significant else ':'

        # Error bars for valid odds ratios
        error = [[odds_ratio - ci.low], [ci.high - odds_ratio]]
        ax.errorbar(
            x=[odds_ratio], y=[y],
            xerr=error, fmt=marker_style,
            color=color, linestyle=line_style,
            markerfacecolor=marker_fill
        )

    # Customize plot
    # ax.axvline(x=1, color='gray', linestyle='--', linewidth=1)
    ax.set_yticks(y_positions)
    ax.set_yticklabels([
        '-LFC (In-Domain vs. Outside-Domain)',
        '+LFC (In-Domain vs. Outside-Domain)',
        '-LFC (High pLDDT vs. Low pLDDT)',
        '+LFC (High pLDDT vs. Low pLDDT)'
    ])
    ax.set_xlabel('Odds Ratio')
    ax.set_title(f'{input_gene} Enrichment Test Odds Ratios')
    plt.tight_layout()
    plt.show()

def enrichment_tests(df_proteinedits, df_coord_struc, workdir, input_gene, pthr, domain_labels):
    # Add the domain column from Uniprot
    df_proteinedits['domain'] = df_coord_struc['domain']
    df_filtered = df_proteinedits[df_proteinedits['mean_Missense_LFC'] != '-'].copy()
    df_filtered['mean_Missense_LFC'] = df_filtered['mean_Missense_LFC'].astype(float)
    df_filtered['mean_Missense_LFC_Z'] = pd.to_numeric(df_filtered['mean_Missense_LFC_Z'], errors='coerce')

    # Split into positive and negative LFC subsets
    positive_df = df_filtered[df_filtered['mean_Missense_LFC'] > 0]
    negative_df = df_filtered[df_filtered['mean_Missense_LFC'] < 0]

    # Results storage
    results = []

    # Perform both tests for positive and negative LFC values
    for df, lfc_type in zip([positive_df, negative_df], ['Positive LFC', 'Negative LFC']):
        # Test 1
        try:
            test1, odds1, ci1 = high_vs_low_pLDDT(data=df, pthr=pthr, domain_labels=domain_labels)
            results.append({
                'lfc_type': lfc_type,
                'test_type': 'High pLDDT vs. Low pLDDT',
                'odds_ratio': test1[0],
                'ci': ci1,
                'p_value': test1[1]
            })
        except (ValueError, AttributeError):
            results.append({
                'lfc_type': lfc_type,
                'test_type': 'High pLDDT vs. Low pLDDT',
                'odds_ratio': np.nan,
                'ci': None,
                'p_value': np.nan
            })

        # Test 2
        try:
            test2, odds2, ci2 = inside_vs_outside_domain(data=df, pthr=pthr)
            results.append({
                'lfc_type': lfc_type,
                'test_type': 'In-Domain vs. Outside-Domain',
                'odds_ratio': test2[0],
                'ci': ci2,
                'p_value': test2[1]
            })
        except (ValueError, AttributeError):
            results.append({
                'lfc_type': lfc_type,
                'test_type': 'In-Domain vs. Outside-Domain',
                'odds_ratio': np.nan,
                'ci': None,
                'p_value': np.nan
            })
    
    plot_enrichment_tests(results=results, pthr=pthr, input_gene=input_gene)