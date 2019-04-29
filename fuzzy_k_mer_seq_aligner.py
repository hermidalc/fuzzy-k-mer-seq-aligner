#!/usr/bin/env python

import numpy as np
from Bio import SeqIO
from argparse import ArgumentParser
from Levenshtein import distance, hamming
from os import path
from pprint import pprint


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

sim_alg_kwargs = {}
if args.sim_alg == 'smith-waterman':
    # Smith-Waterman setup
    sub_mat = substitution_matrix(args.sw_sub_matrix)
    sub_mat_first_chr = list(sub_mat.keys())[0]
    sim_alg_kwargs = {
        'sub_mat': sub_mat,
        'sigma': abs(args.sw_gap_penalty),
        'match_score': sub_mat[sub_mat_first_chr][sub_mat_first_chr],
    }

# parse sequence files to strings
query_seq_rec = SeqIO.read(args.query_seq, args.query_fmt)
target_seq_rec = SeqIO.read(args.target_seq, args.target_fmt)

# load target sequence into fuzzy hashmap
fuzzy_map = {}
k = len(args.fuzzy_seed)
seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '*']
for i in range(len(target_seq_rec) - k + 1):
    target_k_mer_rec = target_seq_rec[i:(i + k)]
    target_k_mer_exact = ''.join([target_k_mer_rec[x]
                                  for x in seed_exact_idxs])
    target_k_mer_fuzzy = ''.join([target_k_mer_rec[x]
                                  for x in seed_fuzzy_idxs])
    if target_k_mer_exact not in fuzzy_map:
        fuzzy_map[target_k_mer_exact] = {target_k_mer_fuzzy: [i]}
    elif target_k_mer_fuzzy not in fuzzy_map[target_k_mer_exact]:
        fuzzy_map[target_k_mer_exact][target_k_mer_fuzzy] = [i]
    else:
        fuzzy_map[target_k_mer_exact][target_k_mer_fuzzy].append(i)

# align query sequence k-mers to target sequence
k_mer_alignments = {}
for j in range(len(query_seq_rec) - k + 1):
    query_k_mer_rec = query_seq_rec[j:(j + k)]
    query_k_mer_exact = ''.join([query_k_mer_rec[x]
                                 for x in seed_exact_idxs])
    query_k_mer_fuzzy = ''.join([query_k_mer_rec[x]
                                 for x in seed_fuzzy_idxs])
    query_k_mer = str(query_k_mer_rec.seq)
    if query_k_mer_exact in fuzzy_map:
        for target_k_mer_fuzzy in fuzzy_map[query_k_mer_exact]:
            if (similarity_score(query_k_mer_fuzzy, target_k_mer_fuzzy,
                                 args.sim_alg, **sim_alg_kwargs) >=
                    args.sim_cutoff):
                target_k_mer_pos = (
                    fuzzy_map[query_k_mer_exact][target_k_mer_fuzzy])
                if query_k_mer in k_mer_alignments:
                    k_mer_alignments[query_k_mer][0].extend(target_k_mer_pos)
                else:
                    k_mer_alignments[query_k_mer] = (target_k_mer_pos, [])
    query_k_mer_rc_rec = query_k_mer_rec.reverse_complement()
    query_k_mer_rc_exact = ''.join([query_k_mer_rc_rec[x]
                                    for x in seed_exact_idxs])
    query_k_mer_rc_fuzzy = ''.join([query_k_mer_rc_rec[x]
                                    for x in seed_fuzzy_idxs])
    query_k_mer_rc = str(query_k_mer_rc_rec.seq)
    if query_k_mer_rc_exact in fuzzy_map:
        for target_k_mer_fuzzy in fuzzy_map[query_k_mer_rc_exact]:
            if (similarity_score(query_k_mer_rc_fuzzy, target_k_mer_fuzzy,
                                 args.sim_alg, **sim_alg_kwargs) >=
                    args.sim_cutoff):
                target_k_mer_rc_pos = (
                    fuzzy_map[query_k_mer_rc_exact][target_k_mer_fuzzy])
                if query_k_mer_rc in k_mer_alignments:
                    k_mer_alignments[query_k_mer_rc][1].extend(
                        target_k_mer_rc_pos)
                else:
                    k_mer_alignments[query_k_mer_rc] = (
                        [], target_k_mer_rc_pos)
