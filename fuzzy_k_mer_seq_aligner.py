#!/usr/bin/env python

from argparse import ArgumentParser
import numpy as np
from operator import itemgetter
from os import path
from pprint import pprint
import pyximport
from Bio import SeqIO
from Bio.Alphabet.IUPAC import (ExtendedIUPACProtein, IUPACAmbiguousDNA,
                                IUPACAmbiguousRNA)
from Bio.Seq import Seq
from Levenshtein import distance, hamming
pyximport.install(language_level=3)
from functions import build_fuzzy_map, iter_kmers


def substitution_matrix(sm_file):
    S = {}
    sm_fh = open(sm_file, 'r')
    s_chs = [''] + sm_fh.readline().strip().split()
    for l in sm_fh.readlines():
        vs = l.strip().split()
        if vs[0] not in S:
            S[vs[0]] = {}
        for i in range(1, len(vs)):
            S[vs[0]][s_chs[i]] = int(vs[i])
    sm_fh.close()
    return S


def smith_waterman_score(m, n, S, sigma):
    F = np.zeros((len(m) + 1, len(n) + 1), dtype=int)
    B = np.zeros((len(m) + 1, len(n) + 1), dtype=int)
    for i in range(1, len(m) + 1):
        for j in range(1, len(n) + 1):
            # END = 0, UP = 1, LEFT = 2, DIAG = 3
            scores = [
                0,
                F[i - 1, j] - sigma,
                F[i, j - 1] - sigma,
                F[i - 1, j - 1] + S[m[i - 1]][n[j - 1]]
            ]
            F[i, j] = max(scores)
            B[i, j] = scores.index(F[i, j])
    i, j = np.unravel_index(np.argmax(F, axis=None), F.shape)
    max_score = str(F[i][j])
    return max_score


def similarity_score(str1, str2, alg, **kwargs):
    if alg == 'levenshtein':
        score = 1 - distance(str1, str2) / len(str1)
    elif alg == 'hamming':
        score = 1 - hamming(str1, str2) / len(str1)
    elif alg == 'smith-waterman':
        score = smith_waterman_score(
            str1, str2, kwargs['sub_mat'], kwargs['sigma']
        ) / (len(str1) * kwargs['match_score'])
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
parser.add_argument('--sim-alg', '-sa', type=str, default='levenshtein',
                    help='string similarity algorithm')
parser.add_argument('--sim-cutoff', '-sc', type=float, default='0.7',
                    help='fuzzy membership similarity cutoff')
parser.add_argument('--sw-sub-matrix', '-sws', type=str,
                    default='NUCLEOTIDE.txt',
                    help='Smith-Waterman substitution matrix file')
parser.add_argument('--sw-gap-penalty', '-swg', type=int, default='1',
                    help='Smith-Waterman gap penalty')
parser.add_argument('--id-cutoff', '-ic', type=float, default='0.7',
                    help='identity cutoff')
args = parser.parse_args()

for file in (args.query_seq, args.target_seq, args.sw_sub_matrix):
    if not path.isfile(file):
        parser.error("File %s doesn't exist" % file)

if args.seq_type == 'dna':
    seq_alphabet = IUPACAmbiguousDNA()
elif args.seq_type == 'rna':
    seq_alphabet = IUPACAmbiguousRNA()
else:
    seq_alphabet = ExtendedIUPACProtein()

sim_alg_kwargs = {}
# Smith-Waterman setup
if args.sim_alg == 'smith-waterman':
    sub_mat = substitution_matrix(args.sw_sub_matrix)
    sub_mat_first_chr = list(sub_mat.keys())[0]
    sim_alg_kwargs = {
        'sub_mat': sub_mat,
        'sigma': abs(args.sw_gap_penalty),
        'match_score': sub_mat[sub_mat_first_chr][sub_mat_first_chr],
    }

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
            if (similarity_score(
                    query_kmer_fuzzy, target_kmer_fuzzy,
                    args.sim_alg, **sim_alg_kwargs) >= args.sim_cutoff):
                for target_kmer_pos in (
                        fuzzy_map[query_kmer_exact][target_kmer_fuzzy]):
                    query_kmer_align_pos['+'].append(
                        (target_kmer_pos, query_kmer_pos))
    if args.seq_type == 'dna':
        query_kmer_rc = str(Seq(query_kmer,
                                alphabet=seq_alphabet).reverse_complement())
        query_kmer_rc_exact = ''.join([query_kmer_rc[x]
                                       for x in seed_exact_idxs])
        query_kmer_rc_fuzzy = ''.join([query_kmer_rc[x]
                                       for x in seed_fuzzy_idxs])
        if query_kmer_rc_exact in fuzzy_map:
            for target_kmer_fuzzy in fuzzy_map[query_kmer_rc_exact]:
                if (similarity_score(
                        query_kmer_rc_fuzzy, target_kmer_fuzzy,
                        args.sim_alg, **sim_alg_kwargs) >= args.sim_cutoff):
                    for target_kmer_pos in (
                            fuzzy_map[query_kmer_rc_exact][target_kmer_fuzzy]):
                        query_kmer_align_pos['-'].append(
                            (target_kmer_pos, query_kmer_pos))
for strand in query_kmer_align_pos:
    # for large lists in-place sort better than sorted
    query_kmer_align_pos[strand].sort(key=itemgetter(0, 1))
