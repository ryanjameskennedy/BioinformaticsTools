#!/usr/bin/env python3

import sys
import argparse

description = '''
------------------------------------------------------------------------------------------------------------------
Searches for ice nucleation motif in fasta or fastq files.
Executed using: python3 ./search_motif.py -f [fasta|fastq] -i <input_filepath> -o <output_filepath>
------------------------------------------------------------------------------------------------------------------
'''

parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-f', help='file format (fasta or fastq)', metavar='FILE_FORMAT', dest='file_format', default='fasta', required=False)
parser.add_argument('--min', help='min amino acid length for intragenic motif space', metavar='MIN_LENGTH', dest='min_length', default=8, type=int,required=False)
parser.add_argument('--max', help='max amino acid length for intragenic motif space', metavar='MAX_LENGTH', dest='max_length', default=20, type=int, required=False)
parser.add_argument('-i', help='input filepath/url', metavar='INPUT', dest='input', nargs='+', required=True)
parser.add_argument('-o', help='output filepath/url', metavar='OUTPUT', dest='output', required=True)
args = parser.parse_args()

import os
import re
import csv
import time
from Bio import SeqIO
from Bio.Seq import Seq

def create_dir(out_fpath):
	folder = "/".join(out_fpath.split("/")[:-1])
	try:
		os.makedirs(folder)
	except OSError:
		print(f"The {folder} has already been created.")

def write_csv(motif_dict, output_filepath):
	with open(output_filepath, 'w+') as csvfile:
		fieldnames = ["DNA sequence", "AA sequence", "Motif count"] #csv header
		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()
		for seq in motif_dict:
			row_dict = {"DNA sequence":seq, "AA sequence":motif_dict[seq][0], "Motif count":motif_dict[seq][1]} #write rows to CSV
			writer.writerow(row_dict)

def translate_seq(seq):
	while len(seq) % 3 != 0:
		seq = seq[:-1] #chop base pair off end to make sequence divisible by 3 (codons)
	return seq.translate()

def search_motif(motif_dict, motif_count_dict, seq):
	trans_seq = translate_seq(seq)
	rep_motif = "GYG"
	seq_split = trans_seq.split(rep_motif)
	if len(seq_split) > 3: #search for repeat motif in translate sequence
		if not any(len(int_space) < args.min_length or len(int_space) > args.max_length for int_space in seq_split[1:-1]):
			if seq in motif_dict:
				motif_dict[seq][1] += 1
			else:
				motif_dict[seq] = [trans_seq, 1] #add to motif dictionary
			for position in range(1, len(trans_seq)-2): #searches the query sequence and stops when there is 26 bp remaining (minimum seed length = 27 bp)
				if trans_seq[position:position+3] == rep_motif:
					dna_motif = seq[(position-1)*3:(position+4)*3]
					aa_motif = trans_seq[position-1:position+4]
					if dna_motif in motif_count_dict:
						motif_count_dict[dna_motif][1] += 1
					else:
						motif_count_dict[dna_motif] = [aa_motif, 1]
	return motif_dict, motif_count_dict

def parse_query(query, file_format):
	motif_dict = {}
	motif_count_dict = {}
	for file in query: #iterate through query files
		with open(file,"r") as handle:
			for record in SeqIO.parse(handle, file_format.lower()): #parses through query sequence(s)
				dna_seq = record.seq #sequence as one line
				dna_seq_rc = dna_seq.reverse_complement() #reverse compliment sequence
				seqs = [dna_seq, dna_seq[1:], dna_seq[2:], dna_seq_rc, dna_seq_rc[1:], dna_seq_rc[2:]] #DNA sequences in 6 frames
				for seq in seqs:
					motif_dict, motif_count_dict = search_motif(motif_dict, motif_count_dict, seq) #translate sequence and search for motif in frames
	motif_count_dict = {seq: seq_count for seq, seq_count in sorted(motif_count_dict.items(), key=lambda item: item[1][1], reverse=True)} #sort dict based on count
	return motif_dict, motif_count_dict

def main():
	create_dir(args.output)
	motif_dict, motif_count_dict = parse_query(args.input, args.file_format)
	write_csv(motif_dict, args.output)
	write_csv(motif_count_dict, args.output.replace(".csv","_motif_counts.csv"))

if __name__ == "__main__":
	t0 = time.time()
	main()
	t1 = time.time()
	print(f"The motif search took {t1-t0} seconds to run.")
