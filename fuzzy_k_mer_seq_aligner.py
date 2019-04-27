#!/usr/bin/env python

import numpy as np
from Bio import SeqIO
from argparse import ArgumentParser
from Levenshtein import distance, hamming
from os import path
from pprint import pprint


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
                    type=str, default='NUCLEOTIDE.txt',
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
query_seq_rec = SeqIO.read(args.query_seq, 'fasta')
ref_seq_rec = SeqIO.read(args.ref_seq, 'fasta')

# load ref sequence into fuzzy hashmap
fuzzy_map = {}
k = len(args.fuzzy_seed)
seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '*']
for i in range(len(ref_seq_rec) - k + 1):
    ref_k_mer_rec = ref_seq_rec[i:(i + k)]
    ref_k_mer_exact = ''.join([ref_k_mer_rec[i] for i in seed_exact_idxs])
    ref_k_mer_fuzzy = ''.join([ref_k_mer_rec[i] for i in seed_fuzzy_idxs])
    if ref_k_mer_exact not in fuzzy_map:
        fuzzy_map[ref_k_mer_exact] = {ref_k_mer_fuzzy: [i]}
    elif ref_k_mer_fuzzy not in fuzzy_map[ref_k_mer_exact]:
        fuzzy_map[ref_k_mer_exact][ref_k_mer_fuzzy] = [i]
    else:
        fuzzy_map[ref_k_mer_exact][ref_k_mer_fuzzy].append(i)

# align query sequence k-mers to ref sequence
k_mer_align_pos = {}
for i in range(len(query_seq_rec) - k + 1):
    query_k_mer_rec = query_seq_rec[i:(i + k)]
    query_k_mer_exact = ''.join([query_k_mer_rec[i] for i in seed_exact_idxs])
    query_k_mer_fuzzy = ''.join([query_k_mer_rec[i] for i in seed_fuzzy_idxs])
    if query_k_mer_exact in fuzzy_map:
        for ref_k_mer_fuzzy in fuzzy_map[query_k_mer_exact]:
            if (similarity_score(query_k_mer_fuzzy, ref_k_mer_fuzzy,
                                 args.sim_alg, **sim_alg_kwargs) >=
                    args.sim_cutoff):
                query_k_mer = str(query_k_mer_rec.seq)
                if query_k_mer in k_mer_align_pos:
                    k_mer_align_pos[query_k_mer].append(i)
                else:
                    k_mer_align_pos[query_k_mer] = [i]
    query_k_mer_rc_rec = query_k_mer_rec.reverse_complement()
    query_k_mer_rc_exact = ''.join(
        [query_k_mer_rc_rec[i] for i in seed_exact_idxs])
    query_k_mer_rc_fuzzy = ''.join(
        [query_k_mer_rc_rec[i] for i in seed_fuzzy_idxs])
    if query_k_mer_rc_exact in fuzzy_map:
        for ref_k_mer_fuzzy in fuzzy_map[query_k_mer_rc_exact]:
            if (similarity_score(query_k_mer_fuzzy, ref_k_mer_fuzzy,
                                 args.sim_alg, **sim_alg_kwargs) >=
                    args.sim_cutoff):
                query_k_mer_rc = str(query_k_mer_rc_rec.seq)
                if query_k_mer_rc in k_mer_align_pos:
                    k_mer_align_pos[query_k_mer_rc].append(i)
                else:
                    k_mer_align_pos[query_k_mer_rc] = [i]
