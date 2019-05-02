#!/usr/bin/env python

from argparse import ArgumentParser
from importlib import import_module
import numpy as np
from operator import itemgetter
from os import path
from pprint import pprint
import pyximport
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Alphabet.IUPAC import (ExtendedIUPACProtein, IUPACAmbiguousDNA,
                                IUPACAmbiguousRNA)
from Bio.Seq import Seq
from Levenshtein import distance, hamming
pyximport.install(language_level=3)
from functions import build_fuzzy_map, iter_kmers


def similarity(query, target, seq_type, algo, aligner):
    if algo == 'levenshtein':
        score = 1 - distance(query, target) / len(target)
    elif algo == 'hamming':
        score = 1 - hamming(query, target) / len(target)
    elif algo == 'smith-waterman':
        if seq_type in ('dna', 'rna'):
            perfect_score = len(target) * aligner.match
        else:
            perfect_score = np.sum([aligner.substitution_matrix[(c, c)]
                                    for c in target])
        # for multiple best scoring alignments get first
        alignment = aligner.align(query, target)[0]
        score = alignment.score / perfect_score
    return score


parser = ArgumentParser()
parser.add_argument('--fuzzy-seed', '-fs', type=str, required=True,
                    help='fuzzy k-mer seed pattern')
parser.add_argument('--query-seq', '-qs', type=str, required=True,
                    help='query sequence file')
parser.add_argument('--query-fmt', '-qf', type=str, default='fasta',
                    help='query sequence format (any Biopython SeqIO '
                         'supported format)')
parser.add_argument('--target-seq', '-ts', type=str, required=True,
                    help='target sequence file')
parser.add_argument('--target-fmt', '-tf', type=str, default='fasta',
                    help='target sequence format (any Biopython SeqIO '
                         'supported format)')
parser.add_argument('--seq-type', '-st', type=str, required=True,
                    choices=['dna', 'rna', 'protein'],
                    help='sequence type')
parser.add_argument('--sim-algo', '-sa', type=str, default='levenshtein',
                    choices=['levenshtein', 'hamming', 'smith-waterman'],
                    help='string similarity algorithm')
parser.add_argument('--sim-cutoff', '-sc', type=float, default='0.7',
                    help='fuzzy membership similarity cutoff')
parser.add_argument('--match-score', '-ms', type=int, default='2',
                    help='match score')
parser.add_argument('--mismatch-score', '-ss', type=int, default='-1',
                    help='mismatch score')
parser.add_argument('--gap-score', '-gs', type=int, default='-1',
                    help='gap score')
parser.add_argument('--sub-matrix', '-sm', type=str, default='blosum62',
                    help='substitution matrix (any Biopython MatrixInfo '
                         'matrix name)')
parser.add_argument('--id-cutoff', '-ic', type=float, default='0.7',
                    help='identity cutoff')
parser.add_argument('--align-fmt', '-af', type=str, default='tabular',
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

# Smith-Waterman setup
if args.sim_algo == 'smith-waterman':
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    if args.seq_type in ('dna', 'rna'):
        aligner.match = args.match_score
        aligner.mismatch = args.mismatch_score
        aligner.gap_score = args.gap_score
    else:
        aligner.substitution_matrix = getattr(
            import_module('Bio.SubsMat.MatrixInfo'), args.sub_matrix)
else:
    aligner = None

# parse sequence files to strings
query_seq_rec = SeqIO.read(
    args.query_seq, args.query_fmt, alphabet=seq_alphabet)
target_seq_rec = SeqIO.read(
    args.target_seq, args.target_fmt, alphabet=seq_alphabet)

# load target sequence into fuzzy hashmap
k = len(args.fuzzy_seed)
seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '*']
fuzzy_map = build_fuzzy_map(str(target_seq_rec.seq), k,
                            seed_exact_idxs, seed_fuzzy_idxs)

# align query sequence k-mers to target sequence
query_kmers = np.empty((len(query_seq_rec) - k + 1,),
                       dtype=np.dtype('U' + str(k)))
query_kmer_align_pos = {'+': [], '-': []}
for query_kmer_pos, query_kmer in iter_kmers(str(query_seq_rec.seq), k):
    query_kmers[query_kmer_pos] = query_kmer
    query_kmer_exact = ''.join([query_kmer[x] for x in seed_exact_idxs])
    query_kmer_fuzzy = ''.join([query_kmer[x] for x in seed_fuzzy_idxs])
    if query_kmer_exact in fuzzy_map:
        for target_kmer_fuzzy in fuzzy_map[query_kmer_exact]:
            similarity_score = similarity(
                query_kmer_fuzzy, target_kmer_fuzzy, args.seq_type,
                args.sim_algo, aligner)
            if similarity_score >= args.sim_cutoff:
                for target_kmer_pos in (
                        fuzzy_map[query_kmer_exact][target_kmer_fuzzy]):
                    query_kmer_align_pos['+'].append(
                        (target_kmer_pos, query_kmer_pos))
    # reverse complement alignment
    if args.seq_type == 'dna':
        query_kmer_rc = str(Seq(query_kmer,
                                alphabet=seq_alphabet).reverse_complement())
        query_kmer_rc_exact = ''.join([query_kmer_rc[x]
                                       for x in seed_exact_idxs])
        query_kmer_rc_fuzzy = ''.join([query_kmer_rc[x]
                                       for x in seed_fuzzy_idxs])
        if query_kmer_rc_exact in fuzzy_map:
            for target_kmer_fuzzy in fuzzy_map[query_kmer_rc_exact]:
                similarity_score = similarity(
                    query_kmer_rc_fuzzy, target_kmer_fuzzy, args.seq_type,
                    args.sim_algo, aligner)
                if similarity_score >= args.sim_cutoff:
                    for target_kmer_pos in (
                            fuzzy_map[query_kmer_rc_exact][target_kmer_fuzzy]):
                        target_kmer_rc_pos = (len(target_seq_rec)
                                              - target_kmer_pos + k)
                        query_kmer_align_pos['-'].append(
                            (target_kmer_rc_pos, query_kmer_pos))

# process query k-mer alignments
for strand in query_kmer_align_pos:
    query_kmer_align_pos[strand].sort(key=itemgetter(0, 1))
    query_kmer_align_pos_grp = []
    for target_seq_pos, query_seq_pos in query_kmer_align_pos[strand]:
        if query_kmer_align_pos_grp:
            if (query_kmer_align_pos_grp[0][0] + len(query_seq_rec)
                    > target_seq_pos + k):
                query_kmer_align_grp_scores = np.zeros(
                    (len(query_seq_rec),), dtype=bool)
                for grp_target_seq_pos, grp_query_seq_pos in (
                        query_kmer_align_pos_grp):
                    query_kmer_align_grp_scores[grp_query_seq_pos:(
                        grp_query_seq_pos + k)] = 1
                query_kmer_align_grp_score = np.sum(
                    query_kmer_align_grp_scores)
                if query_kmer_align_grp_score >= args.id_cutoff:
                    print(strand, query_kmer_align_pos_grp)
                query_kmer_align_pos_grp = [(target_seq_pos, query_seq_pos)]
            else:
                query_kmer_align_pos_grp.append((target_seq_pos,
                                                 query_seq_pos))
        else:
            query_kmer_align_pos_grp.append((target_seq_pos, query_seq_pos))
