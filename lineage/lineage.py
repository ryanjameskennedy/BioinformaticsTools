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
    This script will parse the user-specified MEGAN lineage file extracting the species, genus and phylum level names and place them into a dictionary. 
    This will be used for replacing the species names in the MEGAN species files with the phylum, genus and species names (tab-separated). 

List of classes and functions:
    change_lines, get_lineage_ranks, get_file_ext, write_out, verify_arg.

List of standard modules:
    os, argparse, sys.

List of "non standard" modules:
    ete3

Procedure:
    1. Parse through file and retrieve lineage information using ete3.
    2. Filter undesired rank groups if specified by user and add to out string.
    3. Write out string.

-----------------------------------------------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------------------------------------------

Retrieves lineage on genus and phylum level for MEGAN column. Ranks must be SEPARATED BY SPACES and LOWERCASE.
Executed using: python3 path/to/lineage.py [-r <Rank(s)>] [--filter <Taxonomic_Rank> <Desired_Group>] -i <Input_Filepath> -o <Output_Filepath>
-----------------------------------------------------------------------------------------------------------------------------------------------

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
    '-m',
    '--megan',
    help='MEGAN input file',
    dest='megan',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-r',
    '--ranks',
    help='desired ranks',
    metavar='TAXONOMIC_RANKS',
    dest='ranks',
    nargs='+',
    default=["all"],
    required=False
    )
parser.add_argument(
    '-f',
    '--filter',
    help='filter for group within rank',
    metavar='FILTER',
    dest='filter',
    nargs='+',
    default=False,
    required=False
    )
parser.add_argument(
    '-i',
    help='input file',
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

import os
from ete3 import NCBITaxa
ncbi = NCBITaxa()

def verify_arg():
    if args.input == args.output:
        parser.error('Error: You cannot save the new file under the same name as the MEGAN species file.')
        sys.exit(1)
    if args.filter:
        if len(args.filter) != 2:
            parser.error('Error: The filter flag must be followed by <Taxonomic_Rank> <Desired_Group> e.g. --filter genus Acinetobacter')
            sys.exit(1)
        filter_group_taxid = ncbi.get_name_translator([args.filter[1]])[args.filter[1]][0]
        args.filter.append(filter_group_taxid)
        

def create_dir(out_fpath):
    folder = "/".join(out_fpath.split("/")[:-1])
    try:
        os.makedirs(folder)
    except OSError:
        print(f"The {folder} folder has already been created.")

def write_out(new_lines, output_filepath):
    with open(output_filepath, 'w+') as fout:
        fout.write(new_lines)

def get_file_ext(input_fpath):
    if input_fpath.endswith("csv"):
        separator = ","
    elif input_fpath.endswith("tsv"):
        separator = "\t"
    elif input_fpath.endswith("megan"):
        separator = "\t"
    return separator

def get_taxid(taxa_name):
    try:
        taxid = ncbi.get_name_translator([taxa_name])[taxa_name][0]
    except KeyError:
        print(f'Error: {taxa_name} does not exist in the NCBI database!')
        genus_taxa_name = taxa_name.split(" ")[0]
        try:
            taxid = ncbi.get_name_translator([genus_taxa_name])[genus_taxa_name][0]
        except KeyError:
            print(f'Error: The genus {genus_taxa_name} was removed as it does not exist in the NCBI database!')
            taxid = False
    return taxid

def filter_lineage(taxid_lineage, inv_rank_id_dict, filter_rank, filter_group, filter_group_taxid):
    filter_flag = False
    if filter_rank in inv_rank_id_dict and inv_rank_id_dict[filter_rank] == filter_group_taxid:
        filter_flag = True
    return filter_flag

def edit_rank_col(name_lineage, inv_rank_id_dict, ranks):
    if ranks[0] == "all":
            rank_col = ";".join(name_lineage)
    elif len(ranks) == 1:
        if ranks[0] in inv_rank_id_dict:
            taxid = inv_rank_id_dict[ranks[0]]
            rank_col = ncbi.get_taxid_translator([taxid])[taxid]
        else:
            rank_col = 'NA'
    else:
        for rank in ranks:
            if rank in inv_rank_id_dict:
                taxid = inv_rank_id_dict[rank]
                rank_col = f'{rank_col}{ncbi.get_taxid_translator([taxid])[taxid]};'
            else:
                rank_col = f'{rank_col}NA;'
    return rank_col
    

def get_lineage_ranks(taxa_name):
    #remove = ["candidatus", "uncultured", "bacterium"]
    taxid = get_taxid(taxa_name)
    rank_col = ""
    filter_flag = False
    if taxid:
        taxid_lineage = ncbi.get_lineage(taxid)
        name_lineage_dict = ncbi.get_taxid_translator(taxid_lineage)
        name_lineage = [name_lineage_dict[taxid] for taxid in taxid_lineage]
        rank_id_dict = ncbi.get_rank(taxid_lineage) #convert taxid list to rank dictionary
        inv_rank_id_dict = {rank: taxid for taxid, rank in rank_id_dict.items()} #invert lineage taxid dictionary
        if args.filter:
            filter_flag = filter_lineage(taxid_lineage, inv_rank_id_dict, *args.filter)
        rank_col = edit_rank_col(name_lineage, inv_rank_id_dict, args.ranks)
    return filter_flag, rank_col

def change_lines(input_fpath):
    separator = get_file_ext(input_fpath)
    with open(input_fpath, 'r') as fin:
        new_lines = fin.readline() #header
        for line in fin:
            line = line.rstrip()
            tax_col = line.split(separator)[0]
            taxa_name = tax_col
            if args.megan:
                taxa_name = tax_col.split('"')[1]
            filter_flag, rank_col = get_lineage_ranks(taxa_name)
            if filter_flag:
                new_line = line.replace(tax_col, rank_col)
                new_lines = new_lines + new_line + "\n"
    return new_lines

def main():
    verify_arg()
    
    new_lines = change_lines(args.input)
    create_dir(args.output)
    write_out(new_lines, args.output)

if __name__ == "__main__":
    main()