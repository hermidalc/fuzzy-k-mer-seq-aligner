#!/usr/bin/env python

from argparse import ArgumentParser
from os import path

parser = ArgumentParser()
parser.add_argument('--fuzzy-seed', '-fs',
                    type=str, required=True,
                    help='fuzzy k-mer seed pattern')
parser.add_argument('--query-seq-file', '-qf',
                    type=str, required=True,
                    help='query sequence FASTA file')
parser.add_argument('--ref-seq-file', '-rf',
                    type=str, required=True,
                    help='reference sequence FASTA file')
parser.add_argument('--str-match-alg', '-sa',
                    type=str, default='levenshtein',
                    help='fuzzy string matching algorithm')
args = parser.parse_args()

# parse FASTA sequence files to strings
if path.isfile(args.ref_seq_file):
    ref_seq = ''
    ref_seq_fh = open(args.ref_seq_file)
    ref_seq_fh.readline()
    for line in ref_seq_fh.readlines():
        ref_seq += line.strip()
    ref_seq_fh.close()
else:
    parser.error("File %s doesn't exist" % args.ref_seq_file)
if path.isfile(args.query_seq_file):
    qry_seq = ''
    qry_seq_fh = open(args.query_seq_file)
    qry_seq_fh.readline()
    for line in qry_seq_fh.readlines():
        qry_seq += line.strip()
    qry_seq_fh.close()
else:
    parser.error("File %s doesn't exist" % args.query_seq_file)

# load ref sequence into fuzzy hashmap
fuzzy_dict = {}
k = len(args.fuzzy_seed)
seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
for i in range(len(ref_seq) - k + 1):
    k_mer = ref_seq[i:(i + k)]
    k_mer_exact = ''.join([k_mer[i] for i in seed_exact_idxs])
    k_mer_fuzzy = ''.join([k_mer[i] for i in seed_fuzzy_idxs])
    if k_mer_exact not in fuzzy_dict:
        fuzzy_dict[k_mer_exact] = {k_mer_fuzzy: [i]}
    elif k_mer_fuzzy not in fuzzy_dict[k_mer_exact]:
        fuzzy_dict[k_mer_exact][k_mer_fuzzy] = [i]
    else:
        fuzzy_dict[k_mer_exact][k_mer_fuzzy].append(i)
