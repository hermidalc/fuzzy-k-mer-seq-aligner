#!/usr/bin/env python

from argparse import ArgumentParser
from os import path
from Levenshtein import distance, hamming
import numpy as np
from pprint import pprint

def read_fasta(file):
    seq_fh = open(file)
    seq_fh.readline()
    seq = ''
    for line in seq_fh.readlines():
        seq += line.strip()
    seq_fh.close()
    return seq


def smith_waterman_score(m, n, sigma, sm_file):
    sigma = abs(sigma)
    S = {}
    sub_fh = open(sm_file, 'r')
    s_chs = [''] + sub_fh.readline().strip().split()
    for l in sub_fh.readlines():
        vs = l.strip().split()
        if vs[0] not in S:
            S[vs[0]] = {}
        for i in range(1, len(vs)):
            S[vs[0]][s_chs[i]] = int(vs[i])
    sub_fh.close()
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
        score = 2 * smith_waterman_score(
            str1, str2, kwargs['sigma'], kwargs['sm_file']
        ) / len(str1)
    return score


parser = ArgumentParser()
parser.add_argument('--fuzzy-seed', '-fs',
                    type=str, required=True,
                    help='fuzzy k-mer seed pattern')
parser.add_argument('--query-seq', '-qs',
                    type=str, required=True,
                    help='query sequence FASTA file')
parser.add_argument('--ref-seq', '-rs',
                    type=str, required=True,
                    help='reference sequence FASTA file')
parser.add_argument('--sim-alg', '-sa',
                    type=str, default='levenshtein',
                    help='string similarity algorithm')
parser.add_argument('--sim-cutoff', '-sc',
                    type=float, default='0.7',
                    help='fuzzy membership similarity cutoff')
parser.add_argument('--sw-sub-matrix', '-sws',
                    type=str, default='DNA.txt',
                    help='Smith-Waterman substitution matrix')
parser.add_argument('--sw-gap-penalty', '-swg',
                    type=int, default='1',
                    help='Smith-Waterman gap penalty')
args = parser.parse_args()

for file in (args.query_seq, args.ref_seq, args.sw_sub_matrix):
    if not path.isfile(file):
        parser.error("File %s doesn't exist" % file)
sim_alg_kwargs = {
    'sigma': args.sw_gap_penalty,
    'sm_file': args.sw_sub_matrix,
}

# parse FASTA sequence files to strings
query_seq = read_fasta(args.query_seq)
ref_seq = read_fasta(args.ref_seq)

# load ref sequence into fuzzy hashmap
fuzzy_map = {}
k = len(args.fuzzy_seed)
seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '*']
for i in range(len(ref_seq) - k + 1):
    k_mer = ref_seq[i:(i + k)]
    k_mer_exact = ''.join([k_mer[i] for i in seed_exact_idxs])
    k_mer_fuzzy = ''.join([k_mer[i] for i in seed_fuzzy_idxs])
    if k_mer_exact not in fuzzy_map:
        fuzzy_map[k_mer_exact] = {k_mer_fuzzy: [i]}
    elif k_mer_fuzzy not in fuzzy_map[k_mer_exact]:
        fuzzy_map[k_mer_exact][k_mer_fuzzy] = [i]
    else:
        fuzzy_map[k_mer_exact][k_mer_fuzzy].append(i)

pprint(fuzzy_map)

# align query sequence k-mers to ref sequence
seq_align_pos = {}
for i in range(len(query_seq) - k + 1):
    q_k_mer = query_seq[i:(i + k)]
    q_k_mer_exact = ''.join([k_mer[i] for i in seed_exact_idxs])
    q_k_mer_fuzzy = ''.join([k_mer[i] for i in seed_fuzzy_idxs])
    if q_k_mer_exact in fuzzy_map:
        for r_k_mer_fuzzy in fuzzy_map[q_k_mer_exact]:
            if similarity_score(q_k_mer_fuzzy, r_k_mer_fuzzy, args.sim_alg,
                                **sim_alg_kwargs) >= args.sim_cutoff:
                if q_k_mer in seq_align_pos:
                    seq_align_pos[q_k_mer].append(i)
                else:
                    seq_align_pos[q_k_mer] = []
