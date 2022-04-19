#!/usr/bin/env python

import os
import csv
import json
from Bio import SeqIO
from Bio.Seq import Seq

def create_dir(out_fpath):
	folder = "/".join(out_fpath.split("/")[:-1])
	try:
		os.makedirs(folder)
	except OSError:
		print(f"The {folder} has already been set up.")

def write_csv(count_dict, output_filepath):
    with open(output_filepath, 'w+') as csvfile:
        fieldnames = ["Species", "Condition", "Count"] #header
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for species_name in count_dict:
            for condition, count in count_dict[species_name].items():
                row_dict = {"Species":species_name, "Condition": condition, "Count":count} #write rows to CSV
                writer.writerow(row_dict)

def extract_intron_sites(gff3_fpath):
    intron_sites_dict = {}
    intron_count = 0
    with open(gff3_fpath, 'r') as fin:
        for line in fin:
            if not line.startswith("#"):
                split_line = line.rstrip().split("\t")
                seq_type = split_line[-1]
                if "intron" in seq_type:
                    intron_count += 1
                    seq = split_line[0]
                    start = int(split_line[3])
                    stop = int(split_line[4])
                    strand = split_line[6]
                    if seq in intron_sites_dict:
                        intron_append([strand, start, stop])
                    else:
                        intron_sites_dict[seq] = [[strand, start, stop]]
                        intron_append = intron_sites_dict[seq].append
    return intron_sites_dict, intron_count

def parse_query(file, intron_sites_dict):
    atg_count = 0
    for record in SeqIO.parse(file, "fasta"):
        description = record.description
        contig = description.split(" ")[0]
        if contig in intron_sites_dict:
            for intron_site in intron_sites_dict[contig]:
                strand = intron_site[0]
                start = int(intron_site[1])
                stop = int(intron_site[2])
                seq = str(record.seq)
                intron = seq[start+1:stop+3]
                if strand == "-":
                    intron = str(Seq(intron).reverse_complement()) #reverse complement sequence
                #print(intron)
                for position in range(0, len(intron)-2):
                    if intron[position:position+3] == "ATG":
                        atg_count += 1
    return atg_count

def main():
    #parse_query("../data/pub/fungi/release-48/aspergillus_niger/Aspergillus_niger.ASM285v2.dna.toplevel.fa", 79, 159)
    intron_count_dict = {}
    with open("../data/intron_filepaths.json", 'r') as fin:
        intron_fpaths = json.load(fin)
        for species_name in intron_fpaths:
            fasta_fpath = intron_fpaths[species_name]["fasta"]
            gff3_fpath = intron_fpaths[species_name]["gff3"]
            intron_sites_dict, intron_count = extract_intron_sites(gff3_fpath)
            atg_count = parse_query(fasta_fpath, intron_sites_dict)
            intron_count_dict[species_name] = {"Intron count": intron_count, "ATG count": atg_count}
    create_dir("../data/stats/intron_counts.csv")
    write_csv(intron_count_dict, "../data/stats/intron_counts.csv")

if __name__ == "__main__":
    main()