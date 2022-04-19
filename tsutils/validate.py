#!/usr/bin/python3

import argparse
import sys

description = '''
-----------------------
Title: validate.py
Date: 2020-02-24
Author(s): Ryan Kennedy
-----------------------

Description:
    The script will validate the results of queryIdentifier.py.
List of functions:
    Dictionary, length, for, if, else, rstrip, split.

List of standard modules:
    Json, operator, argparse, sys.

List of "non standard" modules:
    Bio, Bio.Seq.

Procedure:
    1. Iterate through species dictionary generating a link between the accession number and the respective ranks as a list.
    2. Iterate through the result's CSV files to compare the expected species based on NCBI versus the queryIdenetifier.
    3. Write this out to a text file.

---------------------------------------------------------------
'''

usage = '''
---------------------------------------------------------------
Validates the accuracy results from the seedIdentifier program.
Executed using: python3 ./validate.py -c/-n/-p -i <Results_Folder> -o <Output_Filepath>
---------------------------------------------------------------
'''

parser = argparse.ArgumentParser(
                description=description,
                formatter_class=argparse.RawDescriptionHelpFormatter,
                epilog=usage
                )
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 1.0'
    )
parser.add_argument(
    '-c',
    help='sort based on hit count',
    dest='count',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-n',
    help='sort based on normalised hit count (hit count multiplied by inverse proportion of unique seeds)',
    dest='normalised',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-p',
    help='sort based on proportion of hit count from total unique seed count',
    dest='proportion',
    action='store_true',
    required=False
    )
parser.add_argument(
    '-i',
    help='input results folder',
    metavar='INPUT_FOLDER',
    dest='inputFolder',
    default='./Results/',
    required=False
    )
parser.add_argument(
    '-o',
    help='output file path',
    metavar='OUTPUT_FILEPATH',
    dest='outFilepath',
    required=True
    )

args = parser.parse_args()

import os
import csv
import json
import operator
import numpy as np
from matplotlib import pyplot as plt
from ete3 import NCBITaxa
ncbi = NCBITaxa()

def outputWriter(cor_count, csv_error, miss_count, ke_count, in_count, out_count, incor_count, total_count, out_filepath):
    with open(out_filepath, 'w+') as fout:
        output = f"Number of query files with correct top hit: {cor_count}\n"
        output = output + f"Number of query files with insufficient hits: {csv_error}\n"
        output = output + f"Number of query files not in TaxSeedo database: {miss_count}\n"
        output = output + f"Number of KeyErrors: {ke_count}\n"
        output = output + f"Number of species in unique dict: {in_count}\n"
        output = output + f"Number of not in unique dict: {out_count}\n"
        output = output + f"Number of query files with incorrect top hit: {incor_count}\n"
        output = output + f"Number of query files analysed: {cor_count+incor_count+miss_count}\n"
        output = output + f"Total number of query files: {total_count}\n"
        fout.write(output)

def piePlot(count, names, filename):
    pie_chart = plt.figure()
    ax = pie_chart.add_axes([0,0,1,1])
    ax.axis('equal')
    ax.pie(count, labels = names, autopct = '%1.2f%%')
    #plt.show()
    plt.savefig(f"./Stats/{filename}_PC.png")

def getTaxidNames(taxid, accn):
    taxid_lineage = ncbi.get_lineage(taxid)
    rank_id_dict = ncbi.get_rank(taxid_lineage)
    inv_rank_id_dict = {rank: taxid for taxid, rank in rank_id_dict.items()}
    spe_taxid = inv_rank_id_dict['species']
    species = ncbi.get_taxid_translator([spe_taxid])[spe_taxid]
    try:
        gen_taxid = inv_rank_id_dict['genus']
        genus = ncbi.get_taxid_translator([gen_taxid])[gen_taxid]
    except KeyError:
        #print(f"Genus KeyError accession number: {accn}")
        genus = None
    return genus, species

def createValidatorDict(filepath_dictionary, validator_dictionary):
    valid_dict = {}
    with open(filepath_dictionary, 'r') as fin:
        filepath_dict = json.load(fin)
        no_gen_count = 0
        for filepath in filepath_dict:
            accn = filepath_dict[filepath][1] #access accession number from species.json
            taxid = filepath_dict[filepath][2] #access taxid from species.json
            genus, species_name = getTaxidNames(taxid, accn) #get all rank taxids from strain taxid
            if genus == None:
                no_gen_count += 1
                continue
            full_name = filepath_dict[filepath][0].replace("'","").replace("[","").replace("]","")
            species_fullname = full_name.replace(" ", "_").replace("/", "") #edit full name by replacing spaces and forward slashes
            valid_dict[accn] = [genus, species_name, species_fullname] #create database validator dict
    print(f"Number of RefSeq's without a genus level taxid: {no_gen_count}")
    with open(validator_dictionary, 'w+', encoding='utf-8') as fout:
        json.dump(valid_dict, fout)

def errorChecker(incor_hit_list, uniq_dictionary):
    tot_err_count = 0
    ke_count = 0
    pre_count = 0
    with open(uniq_dictionary, 'r') as fin:
        uniq_dict = json.load(fin)
        inv_uniq_dict = list(np.unique(list(uniq_dict.values())))
        for species in incor_hit_list:
            tot_err_count += 1
            try:
                taxid = ncbi.get_name_translator([species])[species][0]
            except KeyError:
                print(species)
                ke_count += 1
                continue
            if taxid in inv_uniq_dict:
                pre_count += 1
    names = ["Not in NCBI taxonomy", "In unique database", "Not in unique database"]
    counts = [ke_count, pre_count, tot_err_count - (pre_count + ke_count)]
    piePlot(counts, names, "incor_hits")
    with open("./Info/incor_hit_list.txt", 'w+') as fout:
        incor_hit_filepaths = "\n".join(map(str, incor_hit_list))+"\n"
        fout.write(incor_hit_filepaths)
    return ke_count, pre_count, tot_err_count - (pre_count + ke_count)

def columnSorter():
    if args.count:
        sort_col = 1
        sort_name = "count"
    elif args.normalised:
        sort_col = 2
        sort_name = "normalised"
    elif args.proportion:
        sort_col = 3
        sort_name = "proportion"
    return sort_col, sort_name

def hitSorter(file):
    sort_col, sort_name = columnSorter()
    with open(file, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader, None)
        sorted_hits = sorted(reader, key=lambda x:float(x[sort_col]), reverse=True)
    return sorted_hits, sort_col

def idValidator(validator_dictionary, proj_val_dictionary, results_folder, out_filepath):
    miss_count = 0
    cor_count = 0
    incor_count = 0
    total_count = 0
    csv_error = 0
    incor_hit_list = []
    with open(validator_dictionary, 'r') as fin, \
        open(proj_val_dictionary, 'r') as fin2:
        valid_dict = json.load(fin)
        proj_val_dict = json.load(fin2)
        rank_dict = {"gen":0,"spe":1,"str":2} #use this dictionary to indicate ranks
        for filename in os.listdir(results_folder):
            if filename.startswith("hits"):
                rank = filename.split("_")[1]
                rank_number = rank_dict[rank] #number to identify expected hit
                file = os.path.join(results_folder, filename)
                exp_accn = "_".join(filename.split("_")[2:4])
                try:
                    exp_hit = valid_dict[exp_accn][rank_number] #uses the rank number to call expected hit based on rank
                except KeyError:
                    try:
                        exp_hit = proj_val_dict[exp_accn][rank_number] #if expected accession not in database validation dict, check project validation dict
                    except KeyError:
                        #print(f"{exp_accn} is not in the TaxSeedo database.")
                        miss_count += 1
                        continue #move to next file
                sorted_hits, sort_col = hitSorter(file)
                total_count += 1
                try:
                    top_hit = sorted_hits[0]
                    top_hit_name = top_hit[0]
                    #top_hit_count = int(top_hit[sort_col])
                    #second_hit = sorted_hits[1]
                    #second_hit_name = second_hit[0]
                    #second_hit_count = int(second_hit[sort_col])
                    #delta_hit = top_hit_count - second_hit_count
                    if exp_hit == top_hit_name:
                        #print(f"File name: {filename}")
                        #print(f"Expected hit: {exp_hit}")
                        #print(f"Top hit: {top_hit_name}")
                        #print(f"Second hit: {second_hit_name}")
                        #print(f"Top hit count: {top_hit_count}\n")
                        cor_count += 1
                    else:
                        incor_hit_list.append(f'./CDSs/refseq/bacteria/{"_".join(filename.split("_")[2:4])}/{filename.replace("hits_spe_","").replace(".csv","_genomic.fna")}')
                        #print("-------------------------Incorrect-------------------------")
                        #print(f"File name: {filename}")
                        #print(f"Expected hit: {exp_hit}")
                        #print(f"Top hit: {top_hit_name}")
                        #print(f"Top hit count: {top_hit_count}")
                        #print(f"Second hit: {second_hit_name}")
                        #print(f"Delta hit count: {delta_hit}")
                        #print("-----------------------------------------------------------\n")
                        incor_count += 1
                except:
                    csv_error += 1
                    print(f"CSV error: {file}")
    names = ["Correct top hit", "Incorrect top hit", "Insufficient hits"]
    counts = [cor_count, incor_count, csv_error]
    piePlot(counts, names, "hit_val")
    ke_count, in_count, out_count = errorChecker(incor_hit_list, "./Uniq/uniq_spe_seeds.json")
    outputWriter(cor_count, csv_error, miss_count, ke_count, in_count, out_count, incor_count, total_count, out_filepath)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    createValidatorDict("./Info/species.json", "./Info/validator.json")
    idValidator("./Info/validator.json", "./Info/proj_validator.json", args.inputFolder, args.outFilepath)
