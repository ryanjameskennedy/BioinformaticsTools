#!/usr/bin/env python3

import argparse
import sys

description = '''
------------------------
Title: cut_ani_dddh.py
Date: 2021-02-20
Author(s): Ryan Kennedy
------------------------
Description:
    This script combines all GGDC output csv files and add ANI scores to respective pairwise comparison.

List of functions:
    extract_info, merge_ani_dddh, import_csv.

List of standard modules:
    csv, os, argparse, sys.

List of "non standard" modules:
    pandas, numpy.

Procedure:
    1. Iterate through files within  project directory containing ANI and GGDC output files.
    2. Parse through lines of ani output file and create a nested dictionary of the pairwise ani scores.
    3. Combine all GGDC output csv files and add ANI scores to respective pairwise comparison.

-----------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------
Extracts ani and dddh scores from files found in specified folders.
Executed using: python3 ./cut_ani_dddh.py -i <Input_Directory> -o <Output_Filepath>
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
    '-d',
    help='data directory',
    metavar='DATA_DIRECTORY',
    dest='data',
    required=False
    )
parser.add_argument(
    '-i',
    help='input directory',
    metavar='INPUT_DIRECTORY',
    dest='input',
    required=False
    )
parser.add_argument(
    '-o',
    help='output filepath',
    metavar='OUTPUT_FILEPATH',
    dest='output',
    required=False
    )

args = parser.parse_args()

import os
import numpy as np
import pandas as pd

def import_csv(dddh_filepath):
    dddh_df = pd.read_csv(dddh_filepath, header=1)
    dddh_df.columns.values[14] = "G+C difference"
    return dddh_df

def ncbi_check(ani_score, dddh_score, query_fname, ref_fname, query_species, ref_species, prefix_list):
    ncbi_flag = False
    if any(prefix in query_fname for prefix in prefix_list) and any(prefix in ref_fname for prefix in prefix_list):
        if ani_score > 96 or dddh_score > 70:
            if query_species != ref_species:
                ncbi_flag = True
        else:
            if query_species == ref_species:
                ncbi_flag = True
    return ncbi_flag

def extract_species_names(data_directory, prefix_list):
    species_names_dict = {}
    for filename in os.listdir(data_directory):
        if any(prefix in filename for prefix in prefix_list):
            file = os.path.join(data_directory, filename)
            with open(file, 'r') as fin:
                header = fin.readline()
                species_name = " ".join(header.split(" ")[1:3])
                species_names_dict[os.path.splitext(filename)[0]] = species_name
        else:
            isolate = os.path.splitext(filename)[0]
            species_names_dict[isolate] = isolate
    return species_names_dict

def merge_ani_dddh(species_names_dict, ani_dict, dddh_fpath_list, output_filepath, prefix_list):
    frames = [import_csv(fpath) for fpath in dddh_fpath_list]
    df = pd.concat(frames)
    df = df.reset_index(drop=True)
    df = df.drop_duplicates()
    query_col = list(df["Query genome"])
    ref_col = list(df["Reference genome"])
    dddh_col = list(df["DDH.1"])
    ani_score_list, query_list, ref_list, rel_list, drop_rows = [], [], [], [], []
    for row_num in range(0, len(query_col)):
        query_fname = query_col[row_num]
        ref_fname = ref_col[row_num]
        query_species = species_names_dict[query_fname]
        ref_species = species_names_dict[ref_fname]
        ani_score = float(ani_dict[query_fname][ref_fname])
        dddh_score = float(dddh_col[row_num])
        ncbi_flag = ncbi_check(ani_score, dddh_score, query_fname, ref_fname, query_species, ref_species, prefix_list)
        if ncbi_flag:
            drop_rows.append(row_num)
            print(f"The pairwise comparison between {query_fname} (species: {query_species}) and {ref_fname} (species: {ref_species}) was removed as a result of their ANI ({ani_score}) and dDDh ({dddh_score}).")
        ani_score_list.append(ani_score)
        query_list.append(query_species)
        ref_list.append(ref_species)
        if query_fname == ref_fname:
            rel_list.append("Intraspecific")
        elif not any(prefix in query_fname for prefix in prefix_list) or not any(prefix in ref_fname for prefix in prefix_list):
            rel_list.append("Isolate")
        elif query_species == ref_species:
            rel_list.append("Intraspecific")
        elif query_species != ref_species:
            query_genus = query_species.split(" ")[0]
            ref_genus = ref_species.split(" ")[0]
            if query_genus == ref_genus:
                rel_list.append("Interspecific")
            else:
                rel_list.append("Intergenus")
    df["ANI"] = ani_score_list
    df["Query genome"] = query_list
    df["Reference genome"] = ref_list
    df["Relationship"] = rel_list
    df.drop(drop_rows, inplace=True)
    df.to_csv(f"{output_filepath}.csv", index=False)
    
def extract_info(data_directory, input_directory, output_filepath):
    current_root = None
    dddh_fpath_list = []
    ani_dict = {}
    prefix_list = ["GCA", "GCF"]
    species_names_dict = extract_species_names(data_directory, prefix_list)
    base = os.path.basename(output_filepath)
    fname = os.path.splitext(base)[0]
    print(fname)
    for root, dirs, files in os.walk(input_directory):
        for filename in files:
            if current_root != root and current_root != None and current_root != input_directory:
                merge_ani_dddh(species_names_dict, ani_dict, dddh_fpath_list, output_filepath, prefix_list)
                dddh_fpath_list = []
                ani_dict = {}
            current_root = root
            if filename.startswith("ggdc_results_"):
                dddh_filepath = os.path.join(root, filename)
                dddh_fpath_list.append(dddh_filepath)
            elif filename.startswith("ani_results_"):
                ani_filepath = os.path.join(root, filename)
                with open(ani_filepath, 'r') as fin:
                    for ANI_line in fin:
                        ANI_line = ANI_line.rstrip()
                        query = ANI_line.split("\t")[0].split("/")[-1].rstrip(".fna").rstrip(".fasta")
                        reference = ANI_line.split("\t")[1].split("/")[-1].rstrip(".fna").rstrip(".fasta")
                        ani_score = ANI_line.split("\t")[2]
                        if query in ani_dict:
                            ani_dict[query][reference] = ani_score
                        else:
                            ani_dict[query] = {reference: ani_score}
    merge_ani_dddh(species_names_dict, ani_dict, dddh_fpath_list, output_filepath, prefix_list)

def main():
    extract_info(args.data, args.input, args.output.rstrip(".csv"))

if __name__ == "__main__":
	main()