#!/usr/bin/env python3

import argparse
import sys

description = '''
------------------------
Title: file_tool.py
Date: 2020-02-25
Author(s): Ryan Kennedy
------------------------
Description:
    This script can unzip .gz files from a user specified folder or individual files (Local input: folder; Cluster input: files) and save them to a user specified folder (output). Alternatively, it can write paths of files in a specified folder (input) to a single file.

List of functions:
    For, if, else, rstrip, split.

List of standard modules:
    Shutil, os, argparse, sys, gzip.

List of "non standard" modules:
    None

Procedure:
    Unzip
    1. Unzip .gz file(s).
    2. Write out in new output folder.
    Get file paths
    1. Parse through folder and save file paths to single file (output)
    Find missing files
    1. Convert lines of filepaths in input files into lists.
    2. Iterate through main filepaths and find those not found in query paths.
    Find string in files
    1. Parse through lines of files within user specified folder.
    2. Out put names of files with query string hit.
    Generate randomly selected list of files
    1. Randomly select filepaths from dictionary.
    2. Output files to user specified folder.
    Generate filepaths text of poorly documented NCBI sequences
    1. Parse through downloaded NCBI database taxids.
    2. Retrieve the lineage dictionary, convert to rank dictionary.
    3. If genus not reported, add to error list.
    Sort ANI output
    1. Parse through ANI results.
    2. Generate a sorted dictionary based on ANI scores.
    3. Convert filepaths to full species name (incl. strain)
    Collect filepaths
    1. Invert species name dictionary (species.json).
    2. Generate filepaths texts from accession numbers as inputs.


-----------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------
Unzips gz files and saves them to specified folder. Writes file paths to file.
Executed using: python3 ./file_tool.py [-ul/-uc/-p/-m/-f/-r/-e/-s/-c/-t] -i <Input> [<Input>] -o <Output>
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
    version='%(prog)s 1.3'
    )
parser.add_argument(
    '-ul',
    '--ulocal',
    help='local execution of unzipping files',
    dest='local',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-uc',
    '--ucluster',
    help='cluster execution of unzipping files',
    dest='cluster',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-p',
    '--path',
    help='get path of files in folder',
    dest='path',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-m',
    '--missing',
    help='get files not found in other folder',
    dest='miss',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-f',
    '--find',
    help='find line in files',
    dest='find',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-r',
    '--random',
    help='generate randomly selected file from list',
    dest='random',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-e',
    '--error',
    help='get filepaths of poorly documented NCBI sequences',
    dest='error',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-s',
    '--sort',
    help='sort information from user specified files found within input folder',
    dest='sort',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-c',
    '--collect',
    help='collect filepaths of input(s)',
    dest='collect',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-d',
    '--dataframe',
    help='create a dataframe for pairwise ANI results',
    dest='dataframe',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-t',
    '--threshold',
    help='threshold for TaxSeedo results',
    dest='threshold',
    required=False
    )
parser.add_argument(
    '-i',
    help='input file/folder',
    metavar='INPUT',
    dest='input',
    nargs='+',
    required=True
    )
parser.add_argument(
    '-o',
    help='output file/folder',
    metavar='OUTPUT',
    dest='output',
    required=True
    )

args = parser.parse_args()

import os
import csv
import gzip
import json
import random
import shutil
from ete3 import NCBITaxa
ncbi = NCBITaxa()

def createFileList(file):
    file_list = []
    with open(file, 'r') as fin:
        for line in fin:
            line = line.rstrip()
            file_list.append(line)
    return file_list

def thresholdCutoff(input_folder, prop_threshold, output_filepath):
    output = ""
    for filename in os.listdir(input_folder):
        file = os.path.join(input_folder, filename)
        with open(file, 'r') as fin:
            for line in fin:
                line = line.rstrip()
                species_name = line.split(",")[0]
                prop = line.split(",")[3]
                if not prop == "Proportion of unique seeds (%)":
                    if prop >= prop_threshold:
                        output = f"{output}{filename}: {species_name}|{prop}\n"
    with open(output_filepath, 'w+') as fout:
        fout.write(output)

def csvWriter(output_dict, output_filepath, header):
    with open(output_filepath, 'w+') as csvfile:
        fieldnames = [""] + sorted(header) #header
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for species in fieldnames:
            try:
                row_dict = {"":species}
                row_dict.update(dict(sorted(output_dict[species].items())))
            except KeyError:
                continue
            writer.writerow(row_dict)

def createDataframe(index_files, output_folder):
    for index_file in index_files:
        new_filename = index_file.split("/")[-1].replace("txt","csv")
        output_filepath = os.path.join(output_folder, new_filename)
        output_dict = {}
        species_list = []
        with open(index_file, 'r') as fin:
            for index_line in fin:
                index_line = index_line.rstrip()
                species = index_line.split("\t")[0]
                species_filepath = index_line.split("\t")[2]
                species_list.append(species)
                with open(species_filepath, 'r') as fin2:
                    for ANI_line in fin2:
                        ANI_line = ANI_line.rstrip()
                        comp_species = ANI_line.split("\t")[0]
                        ANI_score = ANI_line.split("\t")[1]
                        try:
                            output_dict[species].update({comp_species:ANI_score})
                        except:
                            output_dict[species] = {comp_species:ANI_score}
        csvWriter(output_dict, output_filepath, species_list)



def collectFilepaths(filepath_dictionary, species_names, output_filepath):
    output_filepaths = ""
    with open(filepath_dictionary, 'r') as fin:
        filepath_dict = json.load(fin)
        name_dict = {info[1]: filepath for filepath, info in filepath_dict.items()}
        for name in species_names:
            try:
                output_filepaths = output_filepaths + name_dict[name].replace("cds_from_genomic.fna","genomic.fna") + "\n"
            except KeyError:
                print(f"{name} not in species filepaths. Check the name is spelt correctly with the correct strain.")
                continue
    with open(output_filepath, 'w+') as fout:
        fout.write(output_filepaths)

def sortInfo(input_folder, query_file, output_filepath):
    info_dict = {}
    output_file = ""
    with open("Info/species.json", 'r') as fin:
        name_dict = json.load(fin)
        for root, dirs, files in os.walk(input_folder):
            for filename in files:
                if filename == query_file:
                    file = os.path.join(root, filename)
                    first_line = open(file).readline().rstrip()
                    sort_value = float(first_line.split("\t")[2])
                    filepath = first_line.split("\t")[1].replace("/gpfs1/scratch/ryan",".").replace("_genomic.fna","_cds_from_genomic.fna")
                    #accn = name_dict[filepath][1]
                    species_name = name_dict[filepath][0]
                    info_dict[species_name] = sort_value
    sorted_info_dict = {species_name: sort_value for species_name, sort_value in sorted(info_dict.items(), key=lambda item: item[1], reverse=True)}
    with open(output_filepath, 'w+') as fout:
        for ani_result in sorted_info_dict:
            output_file = f"{output_file}{ani_result}\t{sorted_info_dict[ani_result]}\n"
        fout.write(output_file)

def errorChecker(input_dictionary, output_filepath):
    no_genus_list = []
    with open(input_dictionary, 'r') as fin:
        input_dict = json.load(fin)
        for filepath in input_dict:
            taxid = input_dict[filepath][2] #access taxid from species.json
            taxid_lineage = ncbi.get_lineage(taxid)
            rank_id_dict = ncbi.get_rank(taxid_lineage)
            inv_rank_id_dict = {rank: taxid for taxid, rank in rank_id_dict.items()}
            try:
                gen_taxid = inv_rank_id_dict['genus']
                continue
            except KeyError:
                no_genus_list.append(filepath.replace("cds_from_genomic.fna","genomic.fna"))
        no_genus_filepaths = "\n".join(map(str, no_genus_list))+"\n"
    with open(output_filepath, 'w+') as fout:
        fout.write(no_genus_filepaths)

def randomSimulator(input_dictionary, sim_numb, output_folder):
    new_filename = os.path.join(output_folder, f"simulation_{sim_numb}.txt")
    sample_list = []
    with open(input_dictionary, 'r') as fin:
        input_dict = json.load(fin)
        for filepath in input_dict:
            taxid = input_dict[filepath][2] #access taxid from species.json
            taxid_lineage = ncbi.get_lineage(taxid)
            rank_id_dict = ncbi.get_rank(taxid_lineage)
            inv_rank_id_dict = {rank: taxid for taxid, rank in rank_id_dict.items()}
            try:
                gen_taxid = inv_rank_id_dict['genus']
                sample_list.append(filepath)
            except KeyError:
                continue
        random_filepaths = random.sample(sample_list, 1000)
        output_filepaths = "\n".join(map(str, random_filepaths)).replace("cds_from_genomic.fna","genomic.fna")+"\n"
    with open(new_filename, 'w+') as fout:
        fout.write(output_filepaths)

def lineFinder(input_folder, find_term, output_file):
    error_files = ""
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            file = os.path.join(root, filename)
            if not filename.startswith("."):
                with open(file, 'r') as fin:
                    for line in fin:
                        line = line.rstrip()
                        if line.find(find_term) != -1:
                            error_files = error_files + file + "\n"
                            break
    with open(output_file, 'a+') as fout:
        fout.write(error_files)

def findMissing(main_file, query_file, output_file):
    missing_files = ""
    main_list = createFileList(main_file)
    query_list = createFileList(query_file)
    for main_file in main_list:
        if main_file in query_list:
            missing_files = missing_files + main_file + "\n"
    print(missing_files)
    with open(output_file, 'w+') as fout:
        fout.write(missing_files)

def createFilepath(input_folder, output_file):
    filepaths = ""
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            relative_path = os.path.join(root, filename)
            filepaths = filepaths + relative_path + "\n"
    with open(output_file, 'w+') as fout:
        fout.write(filepaths)

def unzipFiles(input_folder, output_folder):
    for root, dirs, files in os.walk(input_folder):
        for filename in files:
            relative_path = os.path.join(root, filename)
            if filename.endswith(".gz"):
                new_filename = filename.rstrip(".gz")
                new_filepath = os.path.join(output_folder, new_filename)
                with gzip.open(filename, 'rb') as fin: #open and read .gz file
                    with open(new_filepath, 'wb') as fout: #open new file path with .gz stripped
                        shutil.copyfileobj(fin, fout) #write to new file with new file name to specified folder

def clusterUnzipFiles(input_file, output_folder):
    filename = input_file.split("/")[-1]
    new_filename = filename.rstrip(".gz")
    new_filepath = os.path.join(output_folder, new_filename)
    with gzip.open(input_file, 'rb') as fin: #open and read .gz file
        with open(new_filepath, 'wb') as fout: #open new file path with .gz stripped
            shutil.copyfileobj(fin, fout) #write to new file with new file name to specified folder

def functionCaller(input, output):
    if args.local:
        if not os.path.exists(output):
            os.mkdir(output)
        unzipFiles(input[0], output)
    elif args.cluster:
        if not os.path.exists(output):
            os.mkdir(output)
        clusterUnzipFiles(input[0], output)
    elif args.path:
        createFilepath(input[0], output)
    elif args.miss:
        findMissing(input[0], input[1], output)
    elif args.find:
        lineFinder(input[0], input[1], output)
    elif args.random:
        numb_of_sims = int(input[1]) + 1
        for sim_numb in range(1, numb_of_sims):
            randomSimulator(input[0], sim_numb, output)
    elif args.error:
        errorChecker(input[0])
    elif args.sort:
        sortInfo(input[0], input[1], output)
    elif args.collect:
        collectFilepaths("./Info/species.json", input, output)
    elif args.dataframe:
        createDataframe(input, output)
    elif args.threshold:
        thresholdCutoff(input[0], args.threshold, output)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    functionCaller(args.input, args.output)
