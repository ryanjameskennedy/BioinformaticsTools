#!/usr/bin/env python3

import argparse
import sys

description = '''
------------------------
Title: lineage_retrieval.py
Date: 2021-03-04
Author(s): Ryan James Kennedy
------------------------

Description:
    This script will parse the user-specified MEGAN lineage file extracting the species, genus and phylum level names and place them into a dictionary. This will be used for replacing the species names in the MEGAN species files with the phylum, genus and species names (tab-separated). 

List of classes and functions:
    change_lines, create_lineage_dict, fetch_lineage, search_rank_lineage_dict, write_out, verify_arg.

List of standard modules:
    os, argparse, sys.

List of "non standard" modules:
    ete3

Procedure:
    1. 

-----------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------
Retrieves lineage on genus and phylum level for MEGAN column. Ranks must be SEPARATED BY SPACES and LOWERCASE.
Executed using: python3 path/to/lineage_retrieval.py [-r <Ranks>] -l <MEGAN_lineage_Filepath> -i <MEGAN_Species_Filepath> -o <New_MEGAN_Filepath>
-----------------------------------------------------------------------------------------------------------
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=usage
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.0.1'
    )
parser.add_argument(
    '-l',
    help='input MEGAN lineage file',
    metavar='LINEAGE_FILE',
    dest='lineage',
    required=True
    )
parser.add_argument(
    '-r',
    '--ranks',
    help='desired ranks',
    metavar='TAXONOMIC_RANKS',
    dest='ranks',
    nargs='+',
    default=["phylum","genus","species"],
    required=False
    )
parser.add_argument(
    '-i',
    help='input MEGAN species file',
    metavar='INPUT_FILE',
    dest='input',
    required=True
    )
parser.add_argument(
    '-o',
    help='output filepath',
    metavar='OUTPUT_FILEPATH',
    dest='output',
    required=True
    )

args = parser.parse_args()

from ete3 import NCBITaxa
ncbi = NCBITaxa()

def verify_arg():
    if args.input == args.output:
        parser.error('Error: You cannot save the new file under the same name as the MEGAN species file.')
        sys.exit(1)
    if args.lineage == args.output:
        parser.error('Error: You cannot save the new file under the same name as the MEGAN lineage file.')
        sys.exit(1)

def write_out(new_lines, output_filepath):
    with open(output_filepath, 'w+') as fout:
        fout.write(new_lines)

def get_lineage_ranks(ranks, species_name):
    rank_col = ""
    taxid = ncbi.get_name_translator(species_name)
    taxid_lineage = ncbi.get_lineage(taxid)
    rank_id_dict = ncbi.get_rank(taxid_lineage) #convert taxid dictionary to rank dictionary
    inv_rank_id_dict = {rank: taxid for taxid, rank in rank_id_dict.items()} #invert lineage taxid dictionary
    for rank in ranks:
        if rank in inv_rank_id_dict:
            rank_col = f'{rank_col}"{inv_rank_id_dict[rank]}"\t'
        else:
            rank_col = f'{rank_col}"NA"\t'
    return rank_col

def search_rank_lineage_dict(rank, inv_rank_lineage_dict):
    if rank in inv_rank_lineage_dict:
        taxid = inv_rank_lineage_dict[rank]
        rank_name = ncbi.get_taxid_translator([taxid])[taxid]
    else:
        rank_name = None
    return rank_name

def fetch_lineage(lineage_dict, lineage_list):
    taxid_lineage_list = [taxid[0] for taxid in ncbi.get_name_translator(lineage_list).values()]
    rank_lineage_dict = ncbi.get_rank(taxid_lineage_list) #convert taxid dictionary to rank dictionary
    inv_rank_lineage_dict = {rank: taxid for taxid, rank in rank_lineage_dict.items()} #invert lineage taxid dictionary
    spe_rank_name = lineage_list[-2] #search_rank_lineage_dict("species", inv_rank_lineage_dict)
    gen_rank_name = search_rank_lineage_dict("genus", inv_rank_lineage_dict)
    phy_rank_name = search_rank_lineage_dict("phylum", inv_rank_lineage_dict)
    if gen_rank_name == None:
        gen_rank_name = "NA"
    if phy_rank_name == None:
        phy_rank_name = "NA"
    lineage_dict[spe_rank_name] = f'"{phy_rank_name}"\t"{gen_rank_name}"\t"{spe_rank_name}"'
    return lineage_dict

def create_lineage_dict(lineage_fpath):
    lineage_dict = {}
    with open(lineage_fpath, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                line = line.rstrip()
                lineage_list = line.split("\t")[0].split('"')[1].split(";")
                lineage_dict = fetch_lineage(lineage_dict, lineage_list)
    return lineage_dict

def change_lines(ranks, lineage_dict, input_fpath):
    new_lines = ""
    with open(input_fpath, 'r') as fin:
        for line in fin:
            line = line.rstrip()
            if line.startswith("#"):
                rank_header = '\t'.join(ranks)
                header = line.replace("Datasets", rank_header)
                #header = line.replace("Datasets", "Phylum\tGenus\tSpecies")
                new_lines = new_lines + header + "\n"
            else:
                species_col = line.split("\t")[0]
                species_name = species_col.split('"')[1]
                #rank_col = get_lineage_ranks(ranks, species_name)
                #new_line = line.replace(species_col, rank_col)
                #new_lines = new_lines + new_line + "\n"
                if species_name in lineage_dict:
                    rank_col = lineage_dict[species_name]
                    new_line = line.replace(species_col, rank_col)
                    new_lines = new_lines + new_line + "\n"
                else:
                    print(f"Species name ({species_name}) not in lineage dictionary. Exiting now...")
                    sys.exit(1)
    return new_lines

def main():
    verify_arg()
    lineage_dict = create_lineage_dict(args.lineage)
    new_lines = change_lines(args.ranks, lineage_dict, args.input)
    write_out(new_lines, args.output)

if __name__ == "__main__":
    main()