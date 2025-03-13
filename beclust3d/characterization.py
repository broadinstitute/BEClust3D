"""
File: characterization.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2025-02-23
Description: Characterization of hits and enrichment tests

"""

import pandas as pd
import os
import subprocess
import shutil
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as stats
from scipy.stats import fisher_exact
from pathlib import Path

# def add_annotation(df_coord_struc, workdir, input_gene, input_uniprot, structureid):
#     """
#     Description
#         Adds annotations to the structural data
    
#     Params
#         df_coord_struc: pandas DataFrame
#             DataFrame containing structural data
#         workdir: str
#             Working directory
#         input_gene: str
#             Input gene name
#         input_uniprot: str
#             Input UniProt ID
#         structureid: str
#             Structure ID

#     Returns
#         df_merged: pandas DataFrame
#             Structural data DataFrame merged with new added annotations
#     """

#     # Download UniProt data
#     ffile = f"{input_uniprot}.txt"
#     target_path = os.path.join("drive/MyDrive/BEPipeline", input_gene)
#     os.makedirs(os.path.dirname(target_path), exist_ok=True)
#     subprocess.run(["wget", f"https://rest.uniprot.org/uniprotkb/{ffile}"])

#     downloaded_file = f"./{ffile}"  # File in Colab's local environment
#     if os.path.exists(downloaded_file):
#         target_file = os.path.join(target_path, ffile)
#         if not os.path.exists(target_file): # check if the destination file already exists
#             shutil.move(downloaded_file, target_path)
#             print(f"{downloaded_file} successfully moved to {target_path}")
#         else:
#             print(f"{target_file} already exists. Skipping move.")
#     else:
#         raise FileNotFoundError(f"File {downloaded_file} not found after download!")

#     df_META = pd.DataFrame()
#     df_META["unipos"] = df_coord_struc["unipos"]
#     df_META["unires"] = df_coord_struc["unires"]

#     uniprot_source_filename = f"{workdir}/{input_uniprot}.txt"
#     uniprot_source_file = open(uniprot_source_filename, "r")
#     uniprot_source_lines = uniprot_source_file.readlines()
#     uniprot_source_file.close()

#     uniprot_processed_lines = [idx for idx in uniprot_source_lines if idx[0:2] == "FT"]
#     uniprot_processed_lines_joined = "\n".join(uniprot_processed_lines)

#     uniprot_processed_filename = f"{workdir}/{input_uniprot}_processed.uniprot"
#     uniprot_processed_file = open(uniprot_processed_filename, "w")
#     uniprot_processed_file.writelines(uniprot_processed_lines)
#     uniprot_processed_file.close()

#     arr_domain = []
#     for i in range(0,len(df_META)):
#         unipos = df_META.at[i, 'unipos']
#         uniaa = df_META.at[i, 'unires']

#         # Domains
#         ## domain (code: DOMAIN)
#         domain = '-'
#         for j in range(0, len(uniprot_processed_lines)):
#             featureline = uniprot_processed_lines[j].strip()
#             if "DOMAIN" in featureline:
#                 positionline = featureline.strip()
#                 positioninfo = positionline[21:len(positionline)]
#                 positiontokens = positioninfo.split('..')

#                 notesline = uniprot_processed_lines[j+1].strip()
#                 notesinfo = notesline[21:len(notesline)]
#                 notestoken = notesinfo.split('=')
#                 notes = notestoken[1].strip()
#                 notes = notes[1:len(notes)-1]

#                 if len(positiontokens) == 1:
#                     position = int(positiontokens[0])
#                     if unipos == position:
#                         domain = str(positioninfo) + ":" + str(notes)
#                 elif len(positiontokens) == 2:
#                     start_position = int(positiontokens[0])
#                     end_position = int(positiontokens[1])
#                     if (unipos >= start_position) and (unipos <= end_position):
#                         domain = str(positioninfo) + ":" + str(notes)
#                 elif len(positiontokens) == 0:
#                     print("feature position can't be empty!")
#                     break;
#                 elif len(positiontokens) > 2:
#                     print("feature positions can't be more than two!")
#                     break;
#         arr_domain.append(domain)
#     df_META['domain'] = arr_domain

#     df_META.to_csv(f"{workdir}/{input_gene}_{input_uniprot}_Uniprotf.tsv", sep = '\t', index=False)
#     df_merged = pd.merge(df_coord_struc, df_META, on=['unipos', 'unires'], how='left')
#     # Resave coord_struc_features with added domains
#     df_merged.to_csv(f"{workdir}/{structureid}_coord_struc_features.tsv", sep = '\t', index=False)
#     return df_merged

# # Test 1: High pLDDT vs. Low pLDDT
# def high_vs_low_pLDDT(data, pthr, domain_labels):
#     """
#     Description
#         Perform enrichment test using Fisher's exact test for high vs. low pLDDT values
    
#     Params
#         data: pandas DataFrame
#             DataFrame containing structural data
#         pthr: float
#             p-value threshold
#         domain_labels: dict
#             Dictionary containing domain labels
    
#     Returns
#         ftest: tuple
#             Fisher's exact test result
#         odds_ratio: float
#             Odds ratio
#         ci: tuple
#             Confidence interval
#     """
    
#     below = data[data['mean_Missense_LFC_Z'] < pthr]
#     above = data[data['mean_Missense_LFC_Z'] >= pthr]
#     table = np.array([
#         [len(below[below['pLDDT_dis'].isin(domain_labels['in'])]),
#             len(below[below['pLDDT_dis'].isin(domain_labels['out'])])],
#         [len(above[above['pLDDT_dis'].isin(domain_labels['in'])]),
#             len(above[above['pLDDT_dis'].isin(domain_labels['out'])])]
#     ])
#     ftest = fisher_exact(table, alternative='two-sided')
#     odds_ratio = stats.contingency.odds_ratio(table)
#     ci = odds_ratio.confidence_interval(confidence_level=0.95)
#     return ftest, odds_ratio, ci

# # Test 2: In-Domain vs. Outside-Domain
# def inside_vs_outside_domain(data, pthr):
#     """
#     Description
#         Perform enrichment test using Fisher's exact test for in-domain vs. outside-domain values
    
#     Params
#         data: pandas DataFrame
#             DataFrame containing structural data
#         pthr: float
#             p-value threshold
    
#     Returns
#         ftest: tuple
#             Fisher's exact test result
#         odds_ratio: float
#             Odds ratio
#         ci: tuple
#             Confidence interval
#     """

#     below = data[data['mean_Missense_LFC_Z'] < pthr].reset_index(drop=True)
#     above = data[data['mean_Missense_LFC_Z'] >= pthr].reset_index(drop=True)
#     below_in_domain = below[below['domain'] != '-']
#     below_out_domain = below[below['domain'] == '-']
#     above_in_domain = above[above['domain'] != '-']
#     above_out_domain = above[above['domain'] == '-']
#     table = np.array([
#         [len(below_in_domain), len(below_out_domain)],
#         [len(above_in_domain), len(above_out_domain)]
#     ])
#     ftest = fisher_exact(table, alternative='two-sided')
#     odds_ratio = stats.contingency.odds_ratio(table)
#     ci = odds_ratio.confidence_interval(confidence_level=0.95)
#     return ftest, odds_ratio, ci

# def plot_enrichment_tests(results, pthr, input_gene):
#     """
#     Description
#         Plot enrichment test results
    
#     Params
#         results: list
#             List of results from enrichment tests
#         pthr: float
#             p-value threshold
#         input_gene: str
#             Input gene name
#     """

#     fig, ax = plt.subplots(figsize=(10, 6))
#     y_positions = [1, 2, 3, 4]  # Fixed y-axis positions for the four lines

#     for i, result in enumerate(results):
#         odds_ratio = result['odds_ratio']
#         ci = result['ci']
#         y = y_positions[i]
#         color = 'red' if i % 2 == 0 else 'blue'
#         if np.isnan(odds_ratio) or np.isinf(odds_ratio):
#             # Placeholder for NaN or infinite odds ratio
#             x_min, x_max = ax.get_xlim()
#             x_mid = (x_min + x_max) / 2
#             ax.plot([x_mid], [y], 'o', color="grey", markerfacecolor='none', linestyle='None')
#             continue

#         # Determine styling
#         is_significant = result['p_value'] <= pthr
#         marker_style = 'o' if is_significant else 'o'
#         marker_fill = color if is_significant else 'none'
#         line_style = '-' if is_significant else ':'

#         # Error bars for valid odds ratios
#         error = [[odds_ratio - ci.low], [ci.high - odds_ratio]]
#         ax.errorbar(
#             x=[odds_ratio], y=[y],
#             xerr=error, fmt=marker_style,
#             color=color, linestyle=line_style,
#             markerfacecolor=marker_fill
#         )

#     # Customize plot
#     ax.set_yticks(y_positions)
#     ax.set_yticklabels([
#         '-LFC (In-Domain vs. Outside-Domain)',
#         '+LFC (In-Domain vs. Outside-Domain)',
#         '-LFC (High pLDDT vs. Low pLDDT)',
#         '+LFC (High pLDDT vs. Low pLDDT)'
#     ])
#     ax.set_xlabel('Odds Ratio')
#     ax.set_title(f'{input_gene} Enrichment Test Odds Ratios')
#     plt.tight_layout()
#     plt.show()

### XYH FUNCTIONS ###

import scipy.stats as stats
from scipy.stats import fisher_exact

def enrichment_test(dfs, val_cols, pval_cols, feature_col, names, labels, 
                    input_gene, pthr, domain_labels):
    """
    Description
        Perform enrichment tests for structural data and plots the results
    """
    results = []

    # Perform both tests for positive and negative LFC values
    for df, val, pval, name in zip(dfs, val_cols, pval_cols, names):
        try:
            test1, odds1, ci1 = pLDDT_test(data=df, pval_col=pval, feature_col=feature_col, 
                                           pthr=pthr, domain_labels=domain_labels)
            results.append({
                'lfc_type': name,
                'test_type': 'High pLDDT vs. Low pLDDT',
                'odds_ratio': test1[0],
                'ci': ci1,
                'p_value': test1[1]
            })
        except (ValueError, AttributeError):
            results.append({
                'lfc_type': name,
                'test_type': 'High pLDDT vs. Low pLDDT',
                'odds_ratio': np.nan, 'ci': None, 'p_value': np.nan
            })

    plot_enrichment_test(results=results, pthr=pthr, input_gene=input_gene, labels=labels)
    
# Test 1: High pLDDT vs. Low pLDDT
def pLDDT_test(data, pval_col, feature_col, pthr, domain_labels):
    """
    Description
        Perform enrichment test using Fisher's exact test for high vs. low pLDDT values
    """
    
    below = data[data[pval_col].astype(float) < pthr]
    above = data[data[pval_col].astype(float) >= pthr]
    table = np.array([
        [len(below[below[feature_col].isin(domain_labels['in'])]),
            len(below[below[feature_col].isin(domain_labels['out'])])],
        [len(above[above[feature_col].isin(domain_labels['in'])]),
            len(above[above[feature_col].isin(domain_labels['out'])])]
    ])
    ftest = fisher_exact(table, alternative='two-sided')
    odds_ratio = stats.contingency.odds_ratio(table)
    ci = odds_ratio.confidence_interval(confidence_level=0.95)
    return ftest, odds_ratio, ci

def plot_enrichment_test(results, pthr, input_gene, labels):
    """
    Description
        Plot enrichment test results
    """

    fig, ax = plt.subplots(figsize=(8, len(labels)))
    y_positions = [int(i)+1 for i in range(len(labels))]  # Fixed y-axis positions for the four lines

    for i, result in enumerate(results):
        odds_ratio = result['odds_ratio']
        ci = result['ci']
        y = y_positions[i]
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
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Odds Ratio')
    ax.set_title(f'{input_gene} Enrichment Test Odds Ratios')
    plt.tight_layout()
    plt.show()
    # plt.savefig()

import requests
import csv

def fetch_uniprot_domains(uniprot_id, output_file="uniprot_domains.tsv"):
    """
    Fetches domain annotations from UniProt given a Uniprot ID and saves as a TSV file
    """
    # UniProt API URL #
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"

    # FETCH DATA #
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Error fetching data: {response.status_code}")
        return

    data = response.json()
    
    # SEQUENCE #
    sequence = data.get("sequence", {}).get("value", "")
    if not sequence:
        print("No sequence found.")
        return
    
    # DOMAINS #
    features = data.get("features", [])
    domain_map = {}

    for feature in features:
        if feature.get("type") in {"Domain", "Repeat"}:
            domain_name = feature.get("description", "Unknown domain")
            begin = int(feature["location"]["start"]["value"])
            end = int(feature["location"]["end"]["value"])
            
            for pos in range(begin, end + 1):
                domain_map[pos] = domain_name

    # Write to TSV file
    with open(output_file, "w", newline="") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        writer.writerow(["Position", "Residue", "Domain"])  # Header
        
        for i, residue in enumerate(sequence, start=1):
            domain = domain_map.get(i, "None")
            writer.writerow([i, residue, domain])

    print(f"Data saved to {output_file}")
    