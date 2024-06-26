"""
File: calculate_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.3

"""

import pandas as pd
from pathlib import Path

def prioritize_by_conservation(df_str_cons,
        workdir, input_gene, 
        screenid, 
        nRandom=1000, 
): 

    df_str_cons_3daggr = pd.DataFrame()
    df_str_cons_3daggr['unipos'] = df_str_cons['unipos']
    df_str_cons_3daggr['unires'] = df_str_cons['unires']
    
    edits_filedir = Path(workdir + input_gene)
    str_cons_filename = edits_filedir / f"randomized_screendata/{input_gene}_{screenid}_struc_consrv_missenseedits_randomized.tsv"
    df_str_cons_edits = pd.read_csv(str_cons_filename, sep = "\t")
    
    taa_wise_sum_LFC = []
    taa_wise_norm_LFC = []

    for aa in range(0, len(df_str_cons_edits)):
        naa_str = df_str_cons_edits.at[aa,'Naa'] # neighboring amino acids
        naa_pos_str = df_str_cons_edits.at[aa,'Naa_pos'] # neighboring residue positions
        taa_LFC = df_str_cons_edits.at[aa,'mean_missense_LFC'] # target LFC 
        naa_list = naa_str.split(';')
        naa_pos_list = naa_pos_str.split(';')

        taa_naa_wBE_LFC = 0; # residues that are conserved and with a BE edit value
        sum_taa_naa_LFC = 0.0

        # if taa_pos_conservation == "conserved":
        if taa_LFC != '-':
            taa_naa_wBE_LFC = 1
            sum_taa_naa_LFC = float(taa_LFC)

        for j in range(0, len(naa_list)):
            naa_pos = int(naa_pos_list[j])
            # if naa_pos_conservation == "conserved":
            naa_LFC = df_str_cons_edits.at[naa_pos-1,'mean_missense_LFC']
            if naa_LFC != '-':
                sum_taa_naa_LFC = sum_taa_naa_LFC + float(naa_LFC)
                taa_naa_wBE_LFC = taa_naa_wBE_LFC + 1

        if taa_naa_wBE_LFC == 0:
            taa_wise_sum_LFC.append('-')
            taa_wise_norm_LFC.append('-')
        else:
            taa_wise_sum_LFC.append(str(round(sum_taa_naa_LFC, 3))) 
            taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))
    
    header_LFC = screenid + "_LFC"
    header_LFC3D = screenid + "_LFC3D"
    df_str_cons_3daggr[header_LFC] = df_str_cons_edits['mean_missense_LFC']
    df_str_cons_3daggr[header_LFC3D] = taa_wise_norm_LFC
    
    for r in range(0, nRandom):
        col_head = 'mean_missense_LFCr' + str(r+1)
        taa_wise_sum_LFC = []
        taa_wise_norm_LFC = []

        for aa in range(0, len(df_str_cons_edits)):
            naa_str = df_str_cons_edits.at[aa,'Naa'] # neighboring amino acids
            naa_pos_str = df_str_cons_edits.at[aa,'Naa_pos'] # neighboring residue positions
            taa_LFC = df_str_cons_edits.at[aa,col_head] # target LFC
            naa_list = naa_str.split(';')
            naa_pos_list = naa_pos_str.split(';')

            taa_naa_wBE_LFC = 0; # residues that are conserved and with a BE edit value
            sum_taa_naa_LFC = 0.0

            # if taa_pos_conservation == "conserved":
            if taa_LFC != '-':
                taa_naa_wBE_LFC = 1
                sum_taa_naa_LFC = float(taa_LFC)

            for j in range(0, len(naa_list)):
                naa_pos = int(naa_pos_list[j])
                # if naa_pos_conservation == "conserved":
                naa_LFC = df_str_cons_edits.at[naa_pos-1, col_head]
                if naa_LFC != '-':
                    sum_taa_naa_LFC = sum_taa_naa_LFC + float(naa_LFC)
                    taa_naa_wBE_LFC = taa_naa_wBE_LFC + 1

            if taa_naa_wBE_LFC == 0:
                taa_wise_sum_LFC.append('-')
                taa_wise_norm_LFC.append('-')
            else:
                taa_wise_sum_LFC.append(str(round(sum_taa_naa_LFC, 3))) 
                taa_wise_norm_LFC.append(str(round(sum_taa_naa_LFC/taa_naa_wBE_LFC, 3)))

        header_LFC = screenid + "_LFCr" + str(r+1)
        header_LFC3D = screenid + "_LFC3Dr" + str(r+1)
        df_str_cons_3daggr[header_LFC] = df_str_cons_edits[col_head]
        df_str_cons_3daggr[header_LFC3D] = taa_wise_norm_LFC

    out_filename = edits_filedir / f"LFC3D/{input_gene}_LFC_LFC3D_per_Random_LFC3Dr.tsv"
    df_str_cons_3daggr.to_csv(out_filename, sep = '\t', index=False)
