#!/usr/bin/env python

from argparse import ArgumentParser
from os import path
from importlib import import_module
from Bio.Align import PairwiseAligner
from Bio.Alphabet.IUPAC import (ExtendedIUPACProtein, IUPACAmbiguousDNA,
                                IUPACAmbiguousRNA)
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pyximport
pyximport.install(inplace=True)
from functions import build_alignments, build_fuzzy_map


parser = ArgumentParser()
parser.add_argument('--fuzzy-seed', '-fs', type=str, required=True,
                    help='fuzzy k-mer seed pattern')
parser.add_argument('--query-seq', '-qs', type=str, required=True,
                    help='query FASTA sequence file')
parser.add_argument('--target-seq', '-ts', type=str, required=True,
                    help='target FASTA sequence file')
parser.add_argument('--seq-type', '-st', type=str, required=True,
                    choices=['dna', 'rna', 'protein'], help='sequence type')
parser.add_argument('--sim-algo', '-sa', type=str, default='levenshtein',
                    choices=['levenshtein', 'hamming', 'smith-waterman'],
                    help='string similarity algorithm')
parser.add_argument('--sim-cutoff', '-sc', type=float, default='0.7',
                    help='fuzzy membership similarity cutoff')
parser.add_argument('--match-score', '-ms', type=float, default='2.0',
                    help='match score')
parser.add_argument('--mismatch-score', '-ss', type=float, default='-1.0',
                    help='mismatch score')
parser.add_argument('--open-gap-score', '-og', type=float,
                    help='open gap score')
parser.add_argument('--extend-gap-score', '-eg', type=float,
                    help='extend gap score')
parser.add_argument('--sub-matrix', '-sm', type=str, default='blosum62',
                    help='substitution matrix (any Biopython MatrixInfo '
                         'matrix name)')
parser.add_argument('--max-kmer-gap', '-mg', type=int, default='30',
                    help='maximum gap size when grouping kmers')
parser.add_argument('--expect-thres', '-et', type=float, default='10.0',
                    help='expect value threshold')
parser.add_argument('--align-fmt', '-af', type=str, default='pairwise',
                    choices=['tabular', 'pairwise'],
                    help='alignment output format')
args = parser.parse_args()

for file in (args.query_seq, args.target_seq):
    if not path.isfile(file):
        parser.error("File %s doesn't exist" % file)

if args.seq_type == 'dna':
    seq_alphabet = IUPACAmbiguousDNA()
elif args.seq_type == 'rna':
    seq_alphabet = IUPACAmbiguousRNA()
else:
    seq_alphabet = ExtendedIUPACProtein()

if args.sim_algo == 'smith-waterman':
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    if args.seq_type in ('dna', 'rna'):
        aligner.match = args.match_score
        aligner.mismatch = args.mismatch_score
        if not args.open_gap_score:
            args.open_gap_score = -1.0
        if not args.extend_gap_score:
            args.extend_gap_score = -1.0
    else:
        aligner.substitution_matrix = getattr(
            import_module('Bio.SubsMat.MatrixInfo'), args.sub_matrix)
        if not args.open_gap_score:
            args.open_gap_score = -11.0
        if not args.extend_gap_score:
            args.extend_gap_score = -1.0
    aligner.open_gap_score = args.open_gap_score
    aligner.extend_gap_score = args.extend_gap_score
else:
    aligner = None

args.k = len(args.fuzzy_seed)
args.seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
args.seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '*']
query_seq_fh = open(args.query_seq, 'r')
target_seq_fh = open(args.target_seq, 'r')
for target_seq_title, target_seq in SimpleFastaParser(target_seq_fh):
    fuzzy_map = build_fuzzy_map(target_seq, args)
    for query_seq_title, query_seq in SimpleFastaParser(query_seq_fh):
        build_alignments(fuzzy_map, query_seq, target_seq, aligner, args)
    del fuzzy_map
target_seq_fh.close()
query_seq_fh.close()
