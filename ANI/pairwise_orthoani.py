#!/usr/bin/env python3

import argparse
import sys

description = '''
------------------------
Title: pairwise_orthoani.py
Date: 2021-02-22
Author(s): Ryan Kennedy
------------------------
Description:
    This script can calculate ANI for all vs all scenarios of genomes.

List of functions:
    For, if, else, rstrip, split.

List of standard modules:
    os, argparse, sys.

List of "non standard" modules:
    OrthoANI, multiprocessing, Biopython.

Procedure:
    1. Create pairwise combinations of input batch files.
    2. Calculate ANI values for pairwise comparisons.

-----------------------------------------------------------------------------------------------------------
'''

usage = '''
-----------------------------------------------------------------------------------------------------------
Performs OrthoANI using a specified reference list.
Executed using: python3 ./pairwise_orthoani.py -q <Query> -r <Reference> -o <Output>
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
    version='%(prog)s 0.1'
    )
parser.add_argument(
    '-c', '--cores',
    help='number of cores for parallel processing',
    metavar='NUMBER_OF_CORES',
    dest='n_cores',
    default=1,
    type=int,
    required=False
    )
parser.add_argument(
    '-q',
    help='input query',
    metavar='INPUT_QUERY',
    dest='query',
    required=False
    )
parser.add_argument(
    '--ql',
    help='input query list',
    metavar='INPUT_QUERY_LIST',
    dest='query_list',
    required=False
    )
parser.add_argument(
    '-r',
    help='input reference',
    metavar='INPUT_REFERENCE',
    dest='reference',
    required=False
    )
parser.add_argument(
    '--rl',
    help='input reference list',
    metavar='INPUT_REFERENCE_LIST',
    dest='reference_list',
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

import orthoani
import itertools
from Bio.SeqIO import parse
from multiprocessing import Pool

def pairwise_ani(input_pair):
    ani_output = ""
    query, reference = input_pair
    query_genome = parse(query, "fasta")
    reference_genome = parse(reference, "fasta")
    ani = float(orthoani.orthoani(query_genome, reference_genome))*100
    ani_output = f"{ani_output}{query}\t{reference}\t{ani}"
    return ani_output

def get_genomes(**kwargs):
    input_dict = {}
    for input_type in kwargs:
        genome_fpath = kwargs[input_type]["genome"]
        if genome_fpath:
            genome_list = [genome_fpath]
        else:
            genome_list_fpath = kwargs[input_type]["genome_list"]
            with open(genome_list_fpath, 'r') as fin:
                genome_list = [genome_fpath.rstrip() for genome_fpath in fin]
        input_dict[input_type] = genome_list
    pairwise_list = list(itertools.product(input_dict["query"], input_dict["reference"]))
    return pairwise_list

def main():
    pairwise_list = get_genomes(query={"genome": args.query, "genome_list": args.query_list}, reference={"genome": args.reference, "genome_list": args.reference_list})
    with Pool(args.n_cores) as pool:
        ani_list = pool.map(pairwise_ani, pairwise_list) #Execute the task in parallel producing a list of dictionaries containing
    ani_output = '\n'.join(map(str, ani_list)) 
    with open(args.output, 'w+') as fout:
        fout.write(ani_output)

if __name__ == "__main__":
    main()