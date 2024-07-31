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

aamap = {
        'aa': ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'], 
        'max_asa': [129.0, 274.0, 195.0, 193.0, 167.0, 223.0, 225.0, 104.0, 224.0, 197.0, 
                    201.0, 236.0, 224.0, 240.0, 159.0, 155.0, 172.0, 285.0, 263.0, 174.0], 
        'aa3': ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 
                'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'],
        'aa3cap': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 
                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
}
dssp_dict = {'H':'H', 'G':'H', 'I':'H', 'P':'H', 
             'B':'B', 'E':'B', 'S':'C', 'T':'C', }

def query_uniprot(
        input_uniprot, 
        edits_filedir, out_fasta, 
): 
        """
        Description
                A function to query UniProt for the protein sequence

        Params
                edits_filedir: str, required
                        the Path to the main directory
                input_uniprot: str, required
                        the Uniprot ID for a particular protein
                input_gene: str, required
                        the name of the input human gene

        Returns
                filename of the output protein sequence .fasta                
        """
        
        # fasta to list file #
        ffile = input_uniprot + '.fasta'
        if not os.path.exists(edits_filedir / ffile): 
                _ = wget.download(f'https://rest.uniprot.org/uniprotkb/{ffile}', 
                                out=str(edits_filedir))

        uFasta_file = edits_filedir / ffile
        uFasta = open(uFasta_file, "r")
        header = uFasta.readline() # skip header

        uFasta_list = open(out_fasta, "w")
        uFasta_list.write('unipos\tunires\n')
        j = 0
        for fasta_line in uFasta:
                fasta_line = fasta_line.strip()
                for i in range(len(fasta_line)):
                        uFasta_list.write("%d\t%s\n" % (j+1, fasta_line[i]))
                        j += 1
        uFasta.close()
        uFasta_list.close()

        return uFasta_list

def query_af_and_process(
        input_uniprot, input_gene, 
        edits_filedir, structureid, 
): 
        """
        Description
                A function to query AlphaFold for the protein structure

        Params
                input_uniprot: str, required
                        the Uniprot ID for a particular protein
                input_gene: str, required
                        the name of the input human gene
                edits_filedir: str, required
                        the Path to the main directory
                structureid: str, required
                        the name of the AF and uniprot input

        Returns
                af_processed_file
                        filename of the output processed protein structure .pdb
                coord_pLDDT_file
        """

        # fetch alphafold #
        affile = structureid + '.pdb'
        if not os.path.exists(edits_filedir / affile): 
                _ = wget.download(f'https://alphafold.ebi.ac.uk/files/{affile}', out='')
        af_filename = f"AF_{input_uniprot}.pdb"
        os.rename(affile, edits_filedir / af_filename)

        # Preprocess AF to keep Atoms and DSSP #
        af_source_file = open(edits_filedir / af_filename, "r")
        af_source_lines = af_source_file.readlines()
        af_source_file.close()

        af_processed_lines = [idx for idx in af_source_lines if idx[0:4] == "ATOM"]
        af_processed_lines = ["HEADER\n"] + ["CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"] + af_processed_lines
        af_processed_filename = edits_filedir / f"{structureid}_processed.pdb"
        af_processed_file = open(af_processed_filename, "w")
        af_processed_file.writelines(af_processed_lines)
        af_processed_file.close()

        # Preprocess AF to get Coords and Confidence Vals #
        # open alphafold pdb file and read
        alphafold_filename = edits_filedir / f"{structureid}_processed.pdb"
        alphafold_pdb = PandasPdb()
        alphafold_pdb.read_pdb(str(alphafold_filename))
        atom_df = alphafold_pdb.df['ATOM']
        # open fasta list file, combine fasta and pdb
        fastalist_filename = edits_filedir / f"{input_gene}_{input_uniprot}.tsv"
        fasta_df = pd.read_csv(fastalist_filename, sep = '\t')
        # out data file
        coord_pLDDT_filename = edits_filedir / f"{structureid}_coord.tsv"
        coord_pLDDT_file = open(coord_pLDDT_filename, 'w')
        coord_pLDDT_file.write('\t'.join(['unipos', 'unires', 'x_coord', 'y_coord', 
                                          'z_coord', 'bfactor_pLDDT']) + '\n')

        df_aamap = pd.DataFrame(data=aamap)
        output_data = []
        for i in range(0, len(fasta_df)):
                unipos = fasta_df.at[i, 'unipos']
                uniaa = fasta_df.at[i, 'unires']
                residue_entry = atom_df.loc[atom_df['residue_number'] == int(unipos), ]
                ca_entry = residue_entry.loc[residue_entry['atom_name'] == "CA", ]

                if len(ca_entry) == 0: 
                        x_coord, y_coord, z_coord, b_factor = "-", "-", "-", "-"
                elif len(ca_entry) == 1: 
                        # amino acid in structure (3 letter, all capital)
                        aa_at_ca = ca_entry['residue_name'].iloc[0]
                        # get three letter, all capital code for amino acid at uniprot position
                        uni_res = df_aamap[df_aamap.aa == str(uniaa)].aa3cap
                        if aa_at_ca == uni_res.item(): 
                                x_coord  = round(float(ca_entry['x_coord'].iloc[0]), 3)
                                y_coord  = round(float(ca_entry['y_coord'].iloc[0]), 3)
                                z_coord  = round(float(ca_entry['z_coord'].iloc[0]), 3)
                                b_factor = round(float(ca_entry['b_factor'].iloc[0]), 3)
                        else: 
                                print("residue mismatch", aa_at_ca + ":" + uni_res.item())
                                break
                else: 
                        print("PROBLEM - CHECK")
                        break

                output_data_entry = '\t'.join([str(unipos), str(uniaa), str(x_coord), 
                                               str(y_coord), str(z_coord), str(b_factor)])+'\n'
                output_data.append(output_data_entry)
        
        coord_pLDDT_file.writelines(output_data)
        coord_pLDDT_file.close()

        return af_processed_file, coord_pLDDT_file

def parse_dssp(
        input_uniprot, input_gene, 
        edits_filedir, structureid, 
        
): 
        """
        Description
                A function to parse .dssp file for 

        Params
                edits_filedir: str, required
                        the Path to the main directory
                input_uniprot: str, required
                        the Uniprot ID for a particular protein
                input_gene: str, required
                        the name of the input human gene
                structureid: str, required
                        the name of the AF and uniprot input

        Returns
                dssp_output_filename: 
                        filename of the processed .dssp file
        """

        # Parse DSSP #
        df_aamap = pd.DataFrame(data=aamap)
        # open alphafold dssp file and read, convert to dict
        alphafold_dssp_filename = edits_filedir / f"{structureid}_processed.dssp"
        parser = parseDSSP(alphafold_dssp_filename)
        parser.parse()
        pddict = parser.dictTodataframe()
        # collect info for the relevant chain 
        pddict_ch = pddict.loc[pddict['chain'] == 'A']
        pddict_ch = pddict_ch.fillna('-')
        pddict_ch = pddict_ch.replace(r'^\s*$', '-', regex=True)
        
        # open fasta list file
        fastalist_filename = edits_filedir / f"{input_gene}_{input_uniprot}.tsv"
        fasta_df = pd.read_csv(fastalist_filename, sep = '\t')
        
        # out data file
        dssp_output_filename = edits_filedir / f"{structureid}_dssp_parsed.tsv"
        dssp_output_file = open(dssp_output_filename, 'w')
        dssp_output_file.write('\t'.join(['unipos', 'unires', 'SS9', 'SS3', 'ACC', 'RSA', 
                                          'exposure', 'PHI', 'normPHI', 'PSI', 'normPSI']) + '\n')
        output_data = []

        for i in range(0, len(fasta_df)): 
                unipos = fasta_df.at[i, 'unipos']
                uniaa = fasta_df.at[i, 'unires']
                pddict_ch_entry = pddict_ch.loc[pddict_ch['inscode'] == str(unipos), ]

                if len(pddict_ch_entry) == 0:
                        dssp_SS9, dssp_ASA, dssp_Phi, dssp_Psi, dssp_SS3 = '-', '-', '-', '-', '-'
                elif len(pddict_ch_entry) == 1:
                        dssp_SS9 = pddict_ch_entry['struct'].iloc[0]
                        if dssp_SS9 == "-":
                                dssp_SS9 = "L"
                        dssp_SS9 = dssp_SS9.strip()
                        if dssp_SS9 in dssp_dict.keys(): 
                                dssp_SS3 = dssp_dict[dssp_SS9]
                        else: 
                                dssp_SS3 = "C"

                        dssp_ASA = pddict_ch_entry['acc'].iloc[0]
                        Gly_X_Gly = df_aamap[df_aamap.aa == str(uniaa)].max_asa
                        norm_ASA = round(float(dssp_ASA) / float(Gly_X_Gly), 2)

                        if norm_ASA < 0.05:            exposure = "core"
                        elif 0.05 <= norm_ASA < 0.25:  exposure = "buried"
                        elif 0.25 <= norm_ASA < 0.5:   exposure = "medburied"
                        elif 0.5 <= norm_ASA < 0.75:   exposure = "medexposed"
                        else:                          exposure = "exposed"

                        dssp_Phi = pddict_ch_entry['phi'].iloc[0]
                        norm_Phi = round(float(dssp_Phi) / 180.0, 2)
                        dssp_Psi = pddict_ch_entry['psi'].iloc[0]
                        norm_Psi = round(float(dssp_Psi) / 180.0, 2)
                else:
                        print(pddict_ch_entry)
                        break

                output_data_entry = '\t'.join([str(unipos), str(uniaa), str(dssp_SS9), str(dssp_SS3), 
                                               str(dssp_ASA), str(norm_ASA), str(exposure), str(dssp_Phi), 
                                               str(norm_Phi), str(dssp_Psi), str(norm_Psi)])
                output_data.append(output_data_entry)

        output_data_all = "\n".join(output_data)
        dssp_output_file.writelines(output_data_all)
        dssp_output_file.close()
        
        return dssp_output_filename

def count_aa_within_radius(
        edits_filedir, structureid, 
        radius=6.0, 
): 
        """
        Description
                Count the number of residues within [radius] Angstroms
                of the residue of focus

        Params
                edits_filedir: str, required
                        the Path to the main directory
                structureid: str, required
                        the name of the AF and uniprot input
                radius: float, optional
                        the radius in which to count other residues

        Returns
                df_coord: 
                        a dataframe containing the residues within range of the key residue
        """

        # Counts of amino acids at 6 A distance #
        coord_filename = edits_filedir / f"{structureid}_coord.tsv"
        df_coord = pd.read_csv(coord_filename, sep = "\t")

        taa_wise_contact_count, taa_wise_naa, taa_wise_naa_positions = [], [], []
        for taa in range(0, len(df_coord)):
                t_xcoord = df_coord.at[taa,"x_coord"]
                t_ycoord = df_coord.at[taa,"y_coord"]
                t_zcoord = df_coord.at[taa,"z_coord"]
                taa_all_pairwise_dis, taa_all_naas, taa_all_naas_positions = [], [], []
                
                for naa in range(0, len(df_coord)):
                        if taa != naa:
                                n_xcoord = df_coord.at[naa,"x_coord"]
                                n_ycoord = df_coord.at[naa,"y_coord"]
                                n_zcoord = df_coord.at[naa,"z_coord"]

                                dis_x_sqr = (n_xcoord - t_xcoord)**2
                                dis_y_sqr = (n_ycoord - t_ycoord)**2
                                dis_z_sqr = (n_zcoord - t_zcoord)**2                                
                                pairwise_dist = math.sqrt(dis_x_sqr + dis_y_sqr + dis_z_sqr)
                                if pairwise_dist <= radius: 
                                        taa_all_pairwise_dis.append(pairwise_dist)
                                        taa_all_naas_positions.append(str(df_coord.at[naa, "unipos"]))
                                        taa_all_naas.append(df_coord.at[naa, "unires"])
                
                taa_wise_contact_count.append(len(taa_all_pairwise_dis)) 
                taa_wise_naa.append(';'.join(taa_all_naas))
                taa_wise_naa_positions.append(';'.join(taa_all_naas_positions))
        
        df_coord['Naa_count'] = taa_wise_contact_count
        df_coord['Naa'] = taa_wise_naa
        df_coord['Naa_pos'] = taa_wise_naa_positions
        return df_coord

def degree_of_burial_per_res(
        df_dssp, df_coord, 
        edits_filedir, structureid
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

        # Find degree of burial per residues: maxRSA #
        maxRSA = df_dssp["RSA"].max()
        df_dssp['dBurial'] = round(maxRSA - df_dssp["RSA"], 3)
        df_coord_dssp = pd.merge(df_coord, df_dssp, on=["unipos", "unires"])
        aa_wise_cdBurial = []

        for aa in range(0, len(df_coord_dssp)): 
                taa_dBurial = df_coord_dssp.at[aa, 'dBurial']
                naa_str     = df_coord_dssp.at[aa, 'Naa'] # neighboring amino acids
                naa_pos_str = df_coord_dssp.at[aa, 'Naa_pos'] # neighboring amino acid positions
                naa_list, naa_pos_list = naa_str.split(';'), naa_pos_str.split(';')
                sum_taa_naa_dBurial, norm_sum_taa_naa_dBurial = 0, 0

                for j in range(0, len(naa_list)): 
                        naa_pos = naa_pos_list[j]
                        naa_pBurial = round(df_coord_dssp.at[int(naa_pos)-1, "dBurial"], 2)
                        sum_taa_naa_dBurial = sum_taa_naa_dBurial + naa_pBurial
        
                norm_sum_taa_naa_dBurial = round(sum_taa_naa_dBurial / len(naa_list), 2)
                aa_wise_cdBurial.append(round(norm_sum_taa_naa_dBurial * taa_dBurial, 3))
        df_coord_dssp['normSumdBurial'] = aa_wise_cdBurial

        arr_pLDDT_discrete = []
        for i in range(0, len(df_coord_dssp)):
                pLDDT = df_coord_dssp.at[i, 'bfactor_pLDDT']
                if pLDDT < 50:           pLDDT_discrete = 'very low'
                elif 50 <= pLDDT < 70:   pLDDT_discrete = 'low'
                elif 70 <= pLDDT <= 90:  pLDDT_discrete = 'confident'
                else:                    pLDDT_discrete = 'high'
                arr_pLDDT_discrete.append(pLDDT_discrete)
        df_coord_dssp['pLDDT_dis'] = arr_pLDDT_discrete

        coord_dssp_filename = edits_filedir / f"{structureid}_coord_struc_features.tsv"
        df_coord_dssp.to_csv(coord_dssp_filename, sep="\t", index=False)

        return df_coord_dssp, coord_dssp_filename

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
        
        edits_filedir = Path(f'{workdir}/{input_gene}'.strip('/'))
        if not os.path.exists(edits_filedir):
                os.mkdir(edits_filedir)
        
        out_fasta = edits_filedir / f"{input_gene}_{input_uniprot}.tsv"
        uFasta_list = query_uniprot(input_uniprot, edits_filedir, out_fasta)
        af_processed_file, coord_pLDDT_file = query_af_and_process(input_uniprot, input_gene, edits_filedir, structureid)


        ### query dssp with mkdssp
        dssp_output_filename = parse_dssp(input_uniprot, input_gene, edits_filedir, structureid)

        dssp_filename = edits_filedir / f"{structureid}_dssp_parsed.tsv"
        df_dssp = pd.read_csv(dssp_filename, sep = '\t')
        df_coord = count_aa_within_radius(edits_filedir, structureid)

        df_coord_dssp, coord_dssp_filename = degree_of_burial_per_res(df_dssp, df_coord, edits_filedir, structureid)
        return df_coord_dssp
