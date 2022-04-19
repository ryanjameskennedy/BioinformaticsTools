#!/usr/bin/env python3

import argparse
import sys

description = '''
------------------------
Title: ggdc.py
Date: 2021-04-14
Author(s): Ryan Kennedy
------------------------
Description:
    This script can calculate dDDH values from blast XML out put files.

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
Calculates dDDH from blast output.
Executed using: python3 ./ggdc.py -q <Query> -r <Reference> -o <Output>
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
    nargs='+',
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
    nargs='+',
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

import os
import re
import shlex
import shutil
import itertools
import subprocess

def create_dir(out_fpath):
    folder = "/".join(out_fpath.split("/")[:-1])
    try:
        os.makedirs(folder)
    except OSError:
        shutil.rmtree(folder)
        os.makedirs(folder)
        print(f"The previous {folder} folder was overwritten.")

def get_genomes(**kwargs):
    input_dict = {}
    for input_type in kwargs:
        genome_fpaths = kwargs[input_type]["genome"]
        if genome_fpaths:
            genome_list = genome_fpaths
        else:
            genome_list_fpath = kwargs[input_type]["genome_list"]
            with open(genome_list_fpath, 'r') as fin:
                genome_list = [genome_fpath.rstrip() for genome_fpath in fin]
        input_dict[input_type] = genome_list
    pairwise_list = list(itertools.product(input_dict["query"], input_dict["reference"]))
    return pairwise_list

class Distance():
    def __init__(self):
        self._check_blast()
        self.blast_dir = os.path.join("/".join(args.output.split("/")[:-1]), "blast_dir/")

    @staticmethod
    def _check_blast():
        """Checks blastall and blastn are downloaded"""
        try:
            blastall_proc = subprocess.Popen(['blastall', '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            blastall_stdout, blastall_stderr = blastall_proc.communicate()
            if "command not found" in blastall_stderr:
                print("Please make sure Blast Legacy is installed.")
                sys.exit(0)
            version = re.search(r'version (.+)', stderr)
            return version.group(1)
        except Exception:
            return 'unknown'

    def makeblastdb(self, reference):
        if "GCA" in reference:
            ref_name = "_".join(os.path.basename(reference).split("_")[0:2])
        else:
            ref_name = os.path.basename(reference).replace(".fasta", "")
        blastdb_out_fpath = os.path.join(self.blast_dir, ref_name)
        command_line = f"makeblastdb -in {reference} -dbtype nucl -out {blastdb_out_fpath}"
        args = shlex.split(command_line)
        process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
        stdout, stderr = process.communicate()
        return ref_name, blastdb_out_fpath

    def blastn(self, pairwise_list):
        for pair in pairwise_list:
            query, reference = pair[0], pair[1]
            ref_name, blastdb_out_fpath = self.makeblastdb(reference)
            if "GCA" in query:
                query_name = "_".join(os.path.basename(query).split("_")[0:2])
            else:
                query_name = os.path.basename(query).replace(".fasta", "")
            #out_fpath = f"{"/".join(out_fpath.split('/')[:-1])}/blast_results_{os.path.splitext(out_fpath.split("/")[-1])[0]}.xml"
            blast_out_fpath = os.path.join(self.blast_dir, f"blast_results_{query_name}_vs_{ref_name}.xml")
            command_line = f"blastall -p blastn -i {query} -d {blastdb_out_fpath} -m 7 -a 1 -S 3 -e 10 -F 'm D' -b 100000 -o {blast_out_fpath}"
            args = shlex.split(command_line)
            process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = process.communicate()

    def dddh(self, pairwise_list):
        self.blastn(pairwise_list)
        dddh_output = ""
        for filename in os.listdir(self.blast_dir):
            if filename.endswith(".xml"):
                identity, hsp_length = 0, 0
                query = filename.split("blast_results_")[1].split("_vs_")[0]
                reference = filename.split("_vs_")[1].split(".xml")[0]
                xml_fpath = os.path.join(self.blast_dir, filename)
                with open(xml_fpath, 'r') as fin:
                    for line in fin:
                        if "<Hsp_identity>" in line:
                            ident = float(line.split("<Hsp_identity>")[1].split("<")[0])
                            identity += ident
                        elif "<Hsp_align-len>" in line:
                            hsp_len = float(line.split("<Hsp_align-len>")[1].split("<")[0])
                            hsp_length += hsp_len
                dddh = (1 - (identity/hsp_length))
                dddh_output = f"{dddh_output}{query}\t{reference}\t{dddh}\n"
        return dddh_output

def main():
    #create_dir(self.blast_dir)
    pairwise_list = get_genomes(query={"genome": args.query, "genome_list": args.query_list}, reference={"genome": args.reference, "genome_list": args.reference_list})
    ggdc = Distance()
    dddh_output = ggdc.dddh(pairwise_list)
    with open(args.output, 'w+') as fout:
        fout.write(dddh_output)

if __name__ == "__main__":
    main()