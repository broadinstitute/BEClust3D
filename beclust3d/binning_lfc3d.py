"""
File: binning_lfc3d.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 3.4

"""

import pandas as pd
from pathlib import Path

def binning_lfc3d(df_LFC_LFC3D, 
        screenid, 
        workdir, input_human_gene, 
): 

    df_LFC_LFC3D_dis = pd.DataFrame()
    df_LFC_LFC3D_dis['unipos'] = df_LFC_LFC3D['unipos']
    df_LFC_LFC3D_dis['unires'] = df_LFC_LFC3D['unires']
    
    # add this screen's LFC and LFC3D to the output
    LFC_header = screenid + "_LFC"
    LFC3D_header = screenid + "_LFC3D"
    df_LFC_LFC3D_dis[LFC_header] = df_LFC_LFC3D[LFC_header]
    df_LFC_LFC3D_dis[LFC3D_header] = df_LFC_LFC3D[LFC3D_header]
    
    # generate bin thersholds
    df_LFC_LFC3D_nodash = df_LFC_LFC3D.loc[df_LFC_LFC3D[LFC3D_header] != '-', ]
    df_LFC_LFC3D_nodash = df_LFC_LFC3D_nodash.reset_index(drop=True)
    df_LFC_LFC3D_nodash[LFC3D_header] = df_LFC_LFC3D_nodash[LFC3D_header].astype(float)

    df_LFC3D_neg = df_LFC_LFC3D_nodash.loc[df_LFC_LFC3D_nodash[LFC3D_header] < 0, ]
    df_LFC3D_neg = df_LFC3D_neg.reset_index(drop=True)
    print("length of neg: " + str(len(df_LFC3D_neg)))

    df_LFC3D_pos = df_LFC_LFC3D_nodash.loc[df_LFC_LFC3D_nodash[LFC3D_header] > 0, ]
    df_LFC3D_pos = df_LFC3D_pos.reset_index(drop=True)
    print("length of pos: " + str(len(df_LFC3D_pos)))

    df_LFC3D_neg_stats = df_LFC3D_neg[LFC3D_header].describe()
    print(df_LFC3D_neg_stats)
    df_LFC3D_pos_stats = df_LFC3D_pos[LFC3D_header].describe()
    print(df_LFC3D_pos_stats)

    NEG_10p_v = round(df_LFC3D_neg[LFC3D_header].quantile(0.1), 4) # (bottom 10th percentile)
    print('NEG_10p_v', NEG_10p_v)
    POS_90p_v = round(df_LFC3D_pos[LFC3D_header].quantile(0.9), 4) # (top 90th percentile)
    print('POS_90p_v', POS_90p_v)

    NEG_05p_v = round(df_LFC3D_neg[LFC3D_header].quantile(0.05), 4) # (bottom 5th percentile)
    print('NEG_05p_v', NEG_05p_v)
    POS_95p_v = round(df_LFC3D_pos[LFC3D_header].quantile(0.95), 4) # (top 95th percentile)
    print('POS_95p_v', POS_95p_v)

    ### reformat
    
    # binning and wightning 
    arr_LFC3D_discrete = []
    arr_LFC3D_weight = []

    for i in range(0, len(df_LFC_LFC3D)):
        LFC3D = df_LFC_LFC3D.at[i, LFC3D_header]
        if LFC3D == '-':
            LFC3D_discrete = '-'
            LFC3D_weight = 0.0
        else:
            LFC3Df = round(float(LFC3D), 3)

            ### reformat

            if LFC3Df <= NEG_05p_v:
                LFC3D_discrete = 'NEG_05p'
                LFC3D_weight = -0.95
            elif (LFC3Df <= NEG_10p_v) and (LFC3Df > NEG_05p_v):
                LFC3D_discrete = 'NEG_10p'
                LFC3D_weight = -0.9
            elif (LFC3Df <= df_LFC3D_neg_stats['25%']) and (LFC3Df > NEG_10p_v):
                LFC3D_discrete = 'NEG_25p'
                LFC3D_weight = -0.75
            elif (LFC3Df <= df_LFC3D_neg_stats['50%']) and (LFC3Df > df_LFC3D_neg_stats['25%']):
                LFC3D_discrete = 'NEG_50p'
                LFC3D_weight = -0.5
            elif (LFC3Df <= df_LFC3D_neg_stats['75%']) and (LFC3Df > df_LFC3D_neg_stats['50%']):
                LFC3D_discrete = 'NEG_75p'
                LFC3D_weight = -0.25
            elif (LFC3Df <= df_LFC3D_neg_stats['max']) and (LFC3Df > df_LFC3D_neg_stats['75%']):
                LFC3D_discrete = 'NEG_100p'
                LFC3D_weight = -0.05
            elif (LFC3Df >= df_LFC3D_pos_stats['min']) and (LFC3Df < df_LFC3D_pos_stats['25%']):
                LFC3D_discrete = 'POS_0p'
                LFC3D_weight = 0.05
            elif (LFC3Df >= df_LFC3D_pos_stats['25%']) and (LFC3Df < df_LFC3D_pos_stats['50%']):
                LFC3D_discrete = 'POS_25p'
                LFC3D_weight = 0.25
            elif (LFC3Df >= df_LFC3D_pos_stats['50%']) and (LFC3Df < df_LFC3D_pos_stats['75%']):
                LFC3D_discrete = 'POS_50p'
                LFC3D_weight = 0.50
            elif (LFC3Df >= df_LFC3D_pos_stats['75%']) and (LFC3Df < POS_90p_v):
                LFC3D_discrete = 'POS_75p'
                LFC3D_weight = 0.75
            elif (LFC3Df >= POS_90p_v) and (LFC3Df < POS_95p_v):
                LFC3D_discrete = 'POS_90p'
                LFC3D_weight = 0.90
            elif LFC3Df >= POS_95p_v:
                LFC3D_discrete = 'POS_95p'
                LFC3D_weight = 0.95

        arr_LFC3D_discrete.append(LFC3D_discrete)
        arr_LFC3D_weight.append(LFC3D_weight)

    LFC3D_dis_header = screenid + "_LFC3D_dis"
    LFC3D_wght_header = screenid + "_LFC3D_wght"
    
    df_LFC_LFC3D_dis[LFC3D_dis_header] = arr_LFC3D_discrete
    df_LFC_LFC3D_dis[LFC3D_wght_header] = arr_LFC3D_weight

    edits_filedir = Path(workdir + input_human_gene)
    out_filename = edits_filedir / f"LFC3D/{input_human_gene}_LFC_LFC3D_dis_wght_Signal_Only_per_Screen.tsv"
    df_LFC_LFC3D_dis.to_csv(out_filename, sep = '\t', index=False)

    return df_LFC_LFC3D_dis
