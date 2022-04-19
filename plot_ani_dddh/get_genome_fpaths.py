#!/usr/bin/env python3

import argparse
import sys

description = '''
------------------------
Title: get_genome_fpaths.py
Date: 2021-03-03
Author(s): Ryan James Kennedy
------------------------

Description:
    This script will output a text file containing filepaths for all strains pertaining to top hit species.

List of functions:
    get_filepaths, get_top_hits, write_out, create_dir, verify_arg.

List of standard modules:
    Json, csv, os, argparse, sys.

List of "non standard" modules:
    Bio, Bio.Seq.

Procedure:
    1. Parse bacterial genome information file and create filepath dictionary for species linked to filepaths.
    2. Parse taxonomic assignment files and select top two species-level hits.
    3. Retrieve filepaths for all strains pertaining to top hit species from TaxSeedo and 16S rRNA ANI analysis.

--------------------------------------------------------------------------------------------------------------------------------
'''

usage = '''
----------------------------------------------------------------------------------------------------------------------------------
Parse taxonomic assignment files and select top two species-level hits. Retrieves filepaths for all strains pertaining to top hit species.
Executed using: python3 get_genome_fpaths.py -i <16S_ANI_Results_File> <TaxSeedo_Results_File> -o <Output_Filepath>
----------------------------------------------------------------------------------------------------------------------------------
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=usage
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 0.1'
    )
parser.add_argument(
    '-i',
    help='input file (result files from 16S rRNA SILVA ANI and TaxSeedo)',
    metavar='INPUT_FILE',
    dest='input',
    nargs='+',
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

def verify_arg():
    if len(args.input) != 2:
        parser.error('16S rRNA SILVA ANI or TaxSeedo file is missing, add in the following order: -i <16S_ANI_Results_File> <TaxSeedo_Results_File>')
        sys.exit(1)

def create_dir(out_fpath):
    folder = "/".join(out_fpath.split("/")[:-1])
    try:
        os.makedirs(folder)
    except OSError:
        print(f"The {folder} folder has already been created.")

def write_out(fpath_dict, top_hit_list, output_filepath):
    with open(output_filepath, 'w+') as fout:
        for hit in top_hits:
            for species in fpath_dict[hit]:
                fout.write(species)

def parse_results(input_fpaths):
    top_hit_dict = {"taxseedo": [], "rrna_ani": []}
    for results_fpath in input_fpaths:
        filename = results_fpath.split("/")[-1]
        with open(results_fpath, 'r') as fin:
            if filename.startswith("taxseedo_results_"):
                while len(top_hit_dict["taxseedo"]) < 2:
                    for line in fin:
                        line = line.rstrip()
                        species_name = line.split(",")[0]
                        if "Candidatus" not in species_name or "sp." not in species_name
                            if species_name not in top_hit_dict["taxseedo"]:
                                top_hit_dict["taxseedo"].append(species_name)
            else:
                while len(top_hit_dict["rrna_ani"]) < 2:
                    for line in fin:
                        line = line.rstrip()
                        if not line.startswith("#"):
                            species_name = " ".join(line.split("\t")[1].split(";")[-1].split("_")[:2])
                            if "Candidatus" not in species_name or "sp." not in species_name:
                                if species_name not in top_hit_dict["taxseedo"]:
                                    top_hit_dict["rrna_ani"].append(species_name)
    return top_hit_dict

def get_top_hits(input_fpaths):
    top_hit_dict = parse_results(input_fpaths)
    print(top_hit_dict)

def get_filepaths(bacteria_records):
    refseq_dir = "/".join(bacteria_records.split("/")[:-1])
    with open(bacteria_records, 'r') as fin:
        fpath_dict = {}
        for line in fin:
            if not line.startswith("#"):
                line = line.rstrip()
                species_name = line.split("\t")[7]
                #strain = line.split("\t")[8].split("=")[1]
                fpath = os.path.join(refseq_dir, line.split("\t")[-3].split("/")[-1], "_genomic.fna.gz")
                if species_name in fpath_dict:
                    fpath_dict[species_name].append(fpath)
                else:
                    fpath_dict[species_name] = [fpath]
    return fpath_dict

def main():
    verify_arg()
    top_hit_list = get_top_hits(args.input)
    fpath_dict = get_filepaths("/gpfs1/scratch/shruti/RefSeq/assembly_summary.txt")
    fpath_dict = get_filepaths("/gpfs1/db/ncbi/genomes/refseq/feb2020/fasta/bacteria_records.txt")
    write_out(fpath_dict, top_hit_list, args.output)

if __name__ == "__main__":
    main()