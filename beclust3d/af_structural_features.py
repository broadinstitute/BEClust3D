"""
File: af_structural_features.py
Author: Calvin XiaoYang Hu, Surya Kiran Mani, Sumaiya Iqbal
Date: 2024-06-18
Description: Translated from Notebook 1
             
"""

import pandas as pd
from DSSPparser import parseDSSP
from pathlib import Path
import os.path
from biopandas.pdb import PandasPdb
import math
import wget
import warnings

aamap = {
    'A': {'max_asa': 129.0, 'aa3cap': 'ALA'}, 'R': {'max_asa': 247.0, 'aa3cap': 'ARG'}, 
    'N': {'max_asa': 195.0, 'aa3cap': 'ASN'}, 'D': {'max_asa': 193.0, 'aa3cap': 'ASP'}, 
    'C': {'max_asa': 167.0, 'aa3cap': 'CYS'}, 'E': {'max_asa': 223.0, 'aa3cap': 'GLU'}, 
    'Q': {'max_asa': 225.0, 'aa3cap': 'GLN'}, 'G': {'max_asa': 104.0, 'aa3cap': 'GLY'}, 
    'H': {'max_asa': 224.0, 'aa3cap': 'HIS'}, 'I': {'max_asa': 197.0, 'aa3cap': 'ILE'}, 
    'L': {'max_asa': 201.0, 'aa3cap': 'LEU'}, 'K': {'max_asa': 236.0, 'aa3cap': 'LYS'}, 
    'M': {'max_asa': 224.0, 'aa3cap': 'MET'}, 'F': {'max_asa': 240.0, 'aa3cap': 'PHE'}, 
    'P': {'max_asa': 159.0, 'aa3cap': 'PRO'}, 'S': {'max_asa': 155.0, 'aa3cap': 'SER'}, 
    'T': {'max_asa': 172.0, 'aa3cap': 'THR'}, 'W': {'max_asa': 285.0, 'aa3cap': 'TRP'}, 
    'Y': {'max_asa': 263.0, 'aa3cap': 'TYR'}, 'V': {'max_asa': 174.0, 'aa3cap': 'VAL'}, 
}
dssp_dict = {'H':'H', 'G':'H', 'I':'H', 'P':'H',   # alpha-helix, 3-10 helix, pi-helix, polyproline helix
             'B':'B', 'E':'B', 'S':'C', 'T':'C', } # beta-bridge, strand, bend, turn/loop


def query_uniprot(
        edits_filedir, input_uniprot
): 
    """
    Description
        A function to query UniProt for the protein sequence
    Params
        edits_filedir: str, required
            the Path to the main directory
        input_uniprot: str, required
            the Uniprot ID for a particular protein
    Returns
        filename of the output protein sequence .fasta                
    """

    # QUERY FASTA FILE #
    ffile = input_uniprot + '.fasta'
    if not os.path.exists(edits_filedir / ffile): 
        _ = wget.download(f'https://rest.uniprot.org/uniprotkb/{ffile}', 
                          out=str(edits_filedir))

    uFasta_file = edits_filedir / ffile
    return uFasta_file

def parse_uniprot(
    uFasta_file, out_fasta
): 
    """
    Description
        A function to process UniProt .fasta into 
        a list of positions and amino acids.
    Params
        edits_filedir: str, required
            the Path to the main directory
        input_uniprot: str, required
            the Uniprot ID for a particular protein
        input_gene: str, required
            the name of the input human gene
    Returns
        None
    """
    # OPEN INPUT AND OUTPUT FILES #
    uFasta_list = open(out_fasta, "w")
    uFasta_list.write('unipos\tunires\n')

    uFasta = open(uFasta_file, "r")
    header = uFasta.readline() # skip header

    # READ FASTA SEQUENCE, AND WRITES POS, AMINO ACID #
    j = 0
    for fasta_line in uFasta:
        fasta_line = fasta_line.strip()
        for i in range(len(fasta_line)):
            uFasta_list.write("%d\t%s\n" % (j+1, fasta_line[i]))
            j += 1

    uFasta.close()
    uFasta_list.close()
    return None

def query_af(
    edits_filedir, af_filename, structureid
): 
    """
    Description
        A function to query AlphaFold for the protein structure
    Params
        edits_filedir: str, required
            the Path to the main directory
        input_uniprot: str, required
            the Uniprot ID for a particular protein
        structureid: str, required
            the name of the AF and uniprot input
    Returns
        af_filename
            filename of the output protein structure .pdb
    """

    # QUERY ALPHAFOLD #
    affile = structureid + '.pdb'
    if not os.path.exists(edits_filedir / affile): 
        _ = wget.download(f'https://alphafold.ebi.ac.uk/files/{affile}', out=str(edits_filedir))
    os.rename(edits_filedir / affile, edits_filedir / af_filename)
    return af_filename

def parse_af(
    edits_filedir, 
    af_filename, af_processed_filename, 
): 
    """
    Description
        Process AlphaFold structure for all atoms and their information
    Params
        edits_filedir: str, required
            the Path to the main directory
    Returns
        None
    """

    # PREPROCESS AF TO KEEP ATOMS #
    af_file = open(edits_filedir / af_filename, "r")
    af_lines = af_file.readlines()
    af_file.close()

    af_processed_lines = [idx for idx in af_lines if idx[0:4] == "ATOM"]
    af_processed_lines = ["HEADER\n"] + ["CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"] + af_processed_lines ###
    af_processed_file = open(edits_filedir / af_processed_filename, "w")
    af_processed_file.writelines(af_processed_lines)
    af_processed_file.close()

    return None

def parse_coord(
    edits_filedir, 
    af_processed_filename, 
    fastalist_filename, coord_filename, 
): 
    """
    Description
        Take in processed AlphaFold and processed fasta and parse
        coordinates and values
    Params
        edits_filedir: str, required
            the Path to the main directory
    Returns
        None
    """

    # GET COORDS AND CONFID VALS FROM PROCESSED AF #
    alphafold_pdb = PandasPdb()
    alphafold_pdb.read_pdb(str(edits_filedir / af_processed_filename))
    atom_df = alphafold_pdb.df['ATOM']
    fasta_df = pd.read_csv(edits_filedir / fastalist_filename, sep = '\t')
    coord_file = open(edits_filedir / coord_filename, 'w')
    coord_file.write('\t'.join(['unipos', 'unires', 'x_coord', 'y_coord', 'z_coord', 'bfactor_pLDDT']) + '\n')

    # PARSE OUT X Y Z B DATA FROM PROCESSED FASTA, PROCESSED AF #
    output_data = []
    for i in range(0, len(fasta_df)):
        unipos = fasta_df.at[i, 'unipos'] ### not fast, dict better
        uniaa = fasta_df.at[i, 'unires'] ### not fast, dict better
        residue_entry = atom_df.loc[atom_df['residue_number'] == int(unipos), ] ###
        ca_entry = residue_entry.loc[residue_entry['atom_name'] == "CA", ] ###

        if len(ca_entry) == 0: 
            x_coord, y_coord, z_coord, b_factor = "-", "-", "-", "-"
        elif len(ca_entry) == 1: 
            aa_at_ca = ca_entry['residue_name'].iloc[0] ###
            uni_res = aamap[str(uniaa)]['aa3cap']
            if aa_at_ca == uni_res: 
                x_coord  = round(float(ca_entry['x_coord'].iloc[0]), 3)
                y_coord  = round(float(ca_entry['y_coord'].iloc[0]), 3)
                z_coord  = round(float(ca_entry['z_coord'].iloc[0]), 3)
                b_factor = round(float(ca_entry['b_factor'].iloc[0]), 3)
            else: 
                warnings.warn(f"residue mismatch {aa_at_ca}: {uni_res}")
        elif len(ca_entry) > 1: 
            warnings.warn("PROBLEM - CHECK") ###

        output_data_entry = '\t'.join([str(unipos), str(uniaa), str(x_coord), str(y_coord), str(z_coord), str(b_factor)])+'\n'
        output_data.append(output_data_entry)

    del alphafold_pdb, atom_df, fasta_df
    coord_file.writelines(output_data)
    coord_file.close()
    return None

def query_dssp(
                
): 
    return None

def parse_dssp(
        edits_filedir, 
        alphafold_dssp_filename, fastalist_filename, 
        dssp_parsed_filename, 
): 
    """
    Description
        A function to parse .dssp file for burial, phi, psi, etc
    Params
        edits_filedir: str, required
            the Path to the main directory
    Returns
        None
    """

    # PARSE DSSP #
    parser = parseDSSP(edits_filedir / alphafold_dssp_filename)
    parser.parse()
    pddict = parser.dictTodataframe()
    pddict_ch = pddict.loc[pddict['chain'] == 'A']
    pddict_ch = pddict_ch.fillna('-')
    pddict_ch = pddict_ch.replace(r'^\s*$', '-', regex=True)
    
    # READ FASTA AND DSSP, WRITE PROCESSED DSSP #
    fasta_df = pd.read_csv(edits_filedir / fastalist_filename, sep = '\t')
    dssp_output_file = open(edits_filedir / dssp_parsed_filename, 'w')
    dssp_output_file.write('\t'.join(['unipos', 'unires', 'SS9', 'SS3', 'ACC', 'RSA', 
                                      'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI']) + '\n')
    
    output_data = []
    for i in range(len(fasta_df)): 
        unipos = fasta_df.at[i, 'unipos'] ###
        uniaa = fasta_df.at[i, 'unires'] ###
        pddict_ch_entry = pddict_ch.loc[pddict_ch['inscode'] == str(unipos), ] ###

        if len(pddict_ch_entry) == 0:
            dssp_SS9, dssp_ASA, dssp_Phi, dssp_Psi, dssp_SS3 = '-', '-', '-', '-', '-'
        elif len(pddict_ch_entry) == 1:
            dssp_SS9 = pddict_ch_entry['struct'].iloc[0].strip() ###
            if dssp_SS9 == "-":               dssp_SS9 = "L"
            if dssp_SS9 in dssp_dict.keys():  dssp_SS3 = dssp_dict[dssp_SS9]
            else:                             dssp_SS3 = "C"

            dssp_ASA = pddict_ch_entry['acc'].iloc[0]
            Gly_X_Gly = aamap[str(uniaa)]['max_asa']
            norm_ASA = round(float(dssp_ASA) / float(Gly_X_Gly), 2)
            if           norm_ASA < 0.05:  exposure = "core"
            elif 0.05 <= norm_ASA < 0.25:  exposure = "buried"
            elif 0.25 <= norm_ASA < 0.50:  exposure = "medburied"
            elif 0.50 <= norm_ASA < 0.75:  exposure = "medexposed"
            else:                          exposure = "exposed"

            dssp_Phi = pddict_ch_entry['phi'].iloc[0]
            norm_Phi = round(float(dssp_Phi) / 180.0, 2)
            dssp_Psi = pddict_ch_entry['psi'].iloc[0]
            norm_Psi = round(float(dssp_Psi) / 180.0, 2)
        else:
            warnings.warn(pddict_ch_entry)

        out = '\t'.join([str(unipos), str(uniaa), str(dssp_SS9), str(dssp_SS3), str(dssp_ASA), str(norm_ASA), 
                         str(exposure), str(dssp_Phi), str(norm_Phi), str(dssp_Psi), str(norm_Psi)])
        output_data.append(out)

    output_data_all = "\n".join(output_data)
    dssp_output_file.writelines(output_data_all)
    dssp_output_file.close()
    del fasta_df, pddict, pddict_ch
    return None

def count_aa_within_radius(
        edits_filedir, coord_filename, 
        radius=6.0, 
): 
    """
    Description
        Count the number of residues within [radius] Angstroms
        of the focal residue
    Params
        edits_filedir: str, required
            the Path to the main directory
        radius: float, optional
            the radius in which to count other residues
    Returns
        df_coord: 
            a dataframe containing the residues of the protein, 
            and which residues are within [radius] Angstroms
    """

    # COUNT AMINO ACIDS IN 6A DISTANCE AND TEIR IDENTITY #
    df_coord = pd.read_csv(edits_filedir / coord_filename, sep = "\t")

    taa_count, taa_naa, taa_naa_positions = [], [], []
    for taa in range(len(df_coord)): 
        t_xcoord = df_coord.at[taa, "x_coord"]
        t_ycoord = df_coord.at[taa, "y_coord"]
        t_zcoord = df_coord.at[taa, "z_coord"]

        dis_count, naas, naas_positions = 0, [], []
        for naa in range(len(df_coord)): 
            if taa != naa:
                unires = df_coord.at[naa, "unires"]
                unipos = df_coord.at[naa, "unipos"]

                xcoord = df_coord.at[naa, "x_coord"] - t_xcoord
                ycoord = df_coord.at[naa, "y_coord"] - t_ycoord
                zcoord = df_coord.at[naa, "z_coord"] - t_zcoord
                pairwise_dist = math.sqrt((xcoord)**2 + (ycoord)**2 + (zcoord)**2)
                if pairwise_dist <= radius: 
                    dis_count += 1
                    naas.append(unires)
                    naas_positions.append(str(unipos))
        
        taa_count.append(dis_count)
        taa_naa.append(';'.join(naas))
        taa_naa_positions.append(';'.join(naas_positions))
    
    df_coord['Naa_count'] = taa_count
    df_coord['Naa'] = taa_naa
    df_coord['Naa_pos'] = taa_naa_positions
    return df_coord

def degree_of_burial(
        df_dssp, df_coord, 
        edits_filedir, coord_dssp_filename
): 
    """
    Description
        Calculate the degree of burial per residue with maxRSA metric

    Params
        df_dssp: pandas DataFrame, required
            a dataframe containing the processed .dssp information
        df_coord: pandas DataFrame, required
            a dataframe containing the residues within range of the key residue
        edits_filedir: str, required
            the Path to the main directory
        structureid: str, required
            the name of the AF and uniprot input

    Returns
        df_coord_dssp: pandas DataFrame
            a dataframe of coordinates and structural features
        coord_dssp_filename: 
            filename of df_coord_dssp
    """
    maxRSA = df_dssp["RSA"].max()
    df_dssp['dBurial'] = round(maxRSA - df_dssp["RSA"], 3)
    df_coord_dssp = pd.merge(df_coord, df_dssp, on=["unipos", "unires"])

    # CALCULATE DEGREE OF BURIAL PER RESIDUE normSumdBurial AND CATEGORY pLDDT_dis #
    aa_wise_cdBurial = []
    arr_pLDDT_discrete = []
    for i in range(len(df_coord_dssp)): 
        taa_dBurial  = df_coord_dssp.at[i, 'dBurial']
        naa_list     = df_coord_dssp.at[i, 'Naa'].split(';') # neighboring amino acids
        naa_pos_list = df_coord_dssp.at[i, 'Naa_pos'].split(';') # neighboring amino acid positions

        # CALCULATE #
        sum_dBurial = 0
        for naa_pos in naa_pos_list: 
            sum_dBurial += round(df_coord_dssp.at[int(naa_pos)-1, "dBurial"], 2)
        norm_sum_dBurial = round(sum_dBurial / len(naa_list), 2)
        aa_wise_cdBurial.append(round(norm_sum_dBurial * taa_dBurial, 3))

        # CATEGORIZE #
        pLDDT = df_coord_dssp.at[i, 'bfactor_pLDDT']
        if         pLDDT < 50:  pLDDT_discrete = 'very low'
        elif 50 <= pLDDT < 70:  pLDDT_discrete = 'low'
        elif 70 <= pLDDT < 90:  pLDDT_discrete = 'confident'
        else:                   pLDDT_discrete = 'high'
        arr_pLDDT_discrete.append(pLDDT_discrete)

    df_coord_dssp['normSumdBurial'] = aa_wise_cdBurial
    df_coord_dssp['pLDDT_dis'] = arr_pLDDT_discrete
    df_coord_dssp.to_csv(edits_filedir / coord_dssp_filename, sep="\t", index=False)
    return df_coord_dssp

def af_structural_features(
        workdir, 
        input_gene, input_uniprot, structureid, 
): 
    """
    Description
        Queries Uniprot, AlphaFold, and DSSP
        Processes the information for structural features to input into downstream functions

    Params
        workdir: str, required
            the working directory
        input_gene: str, required
            the name of the input human gene
        input_uniprot: str, required
            the Uniprot ID for a particular protein
        structureid: str, required
            the name of the AF and uniprot input

    Returns
        df_coord_dssp: pandas DataFrame
            a dataframe of coordinates and structural features
    """
    
    edits_filedir = Path(workdir + '/' +  input_gene)
    if not os.path.exists(edits_filedir):
        os.mkdir(edits_filedir)
    
    out_fasta = edits_filedir / f"{input_gene}_{input_uniprot}.tsv"
    uFasta_file = query_uniprot(edits_filedir, input_uniprot)
    parse_uniprot(uFasta_file, out_fasta)

    af_filename = f"AF_{input_uniprot}.pdb"
    query_af(edits_filedir, af_filename, structureid)

    fastalist_filename = f"{input_gene}_{input_uniprot}.tsv"
    af_processed_filename = f"{structureid}_processed.pdb"
    coord_filename = f"{structureid}_coord.tsv"
    parse_af(edits_filedir, af_filename, af_processed_filename)
    parse_coord(edits_filedir, af_processed_filename, fastalist_filename, coord_filename)

    ### query dssp with mkdssp ### query_dssp
    alphafold_dssp_filename = f"{structureid}_processed.dssp"
    dssp_parsed_filename = f"{structureid}_dssp_parsed.tsv"
    parse_dssp(edits_filedir, alphafold_dssp_filename, fastalist_filename, dssp_parsed_filename)

    df_dssp = pd.read_csv(edits_filedir / dssp_parsed_filename, sep = '\t')
    df_coord = count_aa_within_radius(edits_filedir, coord_filename)

    coord_dssp_filename = f"{structureid}_coord_struc_features.tsv"
    df_coord_dssp = degree_of_burial(df_dssp, df_coord, edits_filedir, coord_dssp_filename)
    return df_coord_dssp
