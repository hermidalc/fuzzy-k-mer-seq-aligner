# functions.pyx
# cython: language_level=3, boundscheck=False, wraparound=True, cdivision=True

import numpy as np
cimport numpy as np
from math import log
from operator import itemgetter
from pprint import pprint
from Bio.Seq import reverse_complement
from Levenshtein import distance, hamming


# iterate over k-mers in a sequence
def kmers(str seq, unsigned int k):
    cdef unsigned int i
    for i in range(len(seq) - k + 1):
        yield i, seq[i:(i + k)]


# build fuzzy hash map
def build_fuzzy_map(str seq, args):
    cdef dict fuzzy_map = {}
    cdef unsigned int kmer_pos, x
    cdef str kmer, kmer_exact, kmer_fuzzy
    for kmer_pos, kmer in kmers(seq, args.k):
        kmer_exact = ''.join([kmer[x] for x in args.seed_exact_idxs])
        kmer_fuzzy = ''.join([kmer[x] for x in args.seed_fuzzy_idxs])
        if kmer_exact not in fuzzy_map:
            fuzzy_map[kmer_exact] = {kmer_fuzzy: [kmer_pos]}
        elif kmer_fuzzy not in fuzzy_map[kmer_exact]:
            fuzzy_map[kmer_exact][kmer_fuzzy] = [kmer_pos]
        else:
            fuzzy_map[kmer_exact][kmer_fuzzy].append(kmer_pos)
    return fuzzy_map


def calc_similarity_score(str query, str target, str seq_type, str algo,
                          aligner):
    cdef float score, perfect_score
    if algo == 'levenshtein':
        score = 1 - distance(query, target) / len(target)
    elif algo == 'hamming':
        score = 1 - hamming(query, target) / len(target)
    elif algo == 'smith-waterman':
        if seq_type in ('dna', 'rna'):
            perfect_score = len(target) * aligner.match
        else:
            perfect_score = np.sum(
                [aligner.substitution_matrix[(c, c)] for c in target])
        score = aligner.score(query, target) / perfect_score
    else:
        score = 0
    return score


# build local alignment from query k-mer alignment group
def build_alignment(list alignment_grp, str query_seq, str target_seq, aligner,
                    args):
    cdef float raw_score = 0
    cdef float bit_score
    cdef list query_align_strs = []
    cdef list match_align_strs = []
    cdef list target_align_strs = []
    cdef str query_kmer, query_chr, target_kmer, target_chr, target_seq_gap
    cdef str query_alignment, match_alignment, target_alignment
    cdef unsigned int query_seq_pos, query_seq_start_pos
    cdef unsigned int target_seq_pos, target_seq_start_pos
    cdef unsigned int gap_size, g
    cdef int query_seq_curr_pos = -1
    cdef int target_seq_curr_pos = -1
    for target_seq_pos, query_seq_pos in alignment_grp:
        if target_seq_curr_pos >= 0:
            if target_seq_pos < <unsigned int>target_seq_curr_pos:
                target_seq_start_pos = target_seq_curr_pos
            elif target_seq_pos > <unsigned int>target_seq_curr_pos:

                target_seq_start_pos = target_seq_pos
                gap_size = target_seq_pos - target_seq_curr_pos
                target_seq_gap = target_seq[target_seq_pos:(
                    target_seq_pos + gap_size)]
                target_align_strs.append(target_seq_gap)
            else:
                target_seq_start_pos = target_seq_pos
        else:
            target_seq_start_pos = target_seq_pos
        if query_seq_curr_pos >= 0:
            if query_seq_pos < <unsigned int>query_seq_curr_pos:
                query_seq_start_pos = query_seq_curr_pos
            elif query_seq_pos > <unsigned int>query_seq_curr_pos:
                query_seq_start_pos = query_seq_pos
                gap_size = query_seq_pos - query_seq_curr_pos
                raw_score += sum(
                    args.open_gap_score if g == 0 else
                    args.extend_gap_score for g in range(gap_size))
                query_align_strs.append('-' * gap_size)
                match_align_strs.append(' ' * gap_size)
            else:
                query_seq_start_pos = query_seq_pos
        else:
            query_seq_start_pos = query_seq_pos
        target_seq_curr_pos = target_seq_pos + args.k
        query_seq_curr_pos = query_seq_pos + args.k
        target_kmer = target_seq[target_seq_start_pos:target_seq_curr_pos]
        query_kmer = query_seq[query_seq_start_pos:query_seq_curr_pos]
        query_align_strs.append(query_kmer)
        target_align_strs.append(target_kmer)
        for target_chr, query_chr in zip(
                target_kmer, query_kmer):
            if args.seq_type in ('dna', 'rna'):
                if query_chr == target_chr:
                    raw_score += args.match_score
                    match_align_strs.append('|')
                else:
                    raw_score += args.mismatch_score
                    match_align_strs.append(' ')
            else:
                raw_score += aligner.substitution_matrix[
                    (query_chr, target_chr)]
                if query_chr == target_chr:
                    match_align_strs.append('|')
                else:
                    match_align_strs.append(' ')
    bit_score = (args.ka_gapped_l * raw_score - log(args.ka_gapped_k)) / log(2)
    return {'query': ''.join(query_align_strs),
            'match': ''.join(match_align_strs),
            'target': ''.join(target_align_strs),
            'raw_score': raw_score,
            'bit_score': bit_score,
            'e_value': len(query_seq) * len(target_seq) / (2 ** bit_score),
            'p_value': 2 ** (-bit_score)}


def pairwise_align(dict fuzzy_map, str query_seq, str target_seq, aligner,
                   args):
    # align query k-mers to target sequence
    cdef dict query_kmer_alignments = {'+': [], '-': []}
    cdef unsigned int query_kmer_pos, target_kmer_pos
    cdef float similarity_score
    cdef str query_kmer, query_kmer_rc, query_kmer_exact, query_kmer_rc_exact
    cdef str query_kmer_fuzzy, query_kmer_rc_fuzzy, target_kmer_fuzzy
    for query_kmer_pos, query_kmer in kmers(query_seq, args.k):
        query_kmer_exact = ''.join([query_kmer[x]
                                    for x in args.seed_exact_idxs])
        query_kmer_fuzzy = ''.join([query_kmer[x]
                                    for x in args.seed_fuzzy_idxs])
        if query_kmer_exact in fuzzy_map:
            for target_kmer_fuzzy in fuzzy_map[query_kmer_exact]:
                similarity_score = calc_similarity_score(
                    query_kmer_fuzzy, target_kmer_fuzzy, args.seq_type,
                    args.sim_algo, aligner)
                if similarity_score >= args.sim_cutoff:
                    for target_kmer_pos in (fuzzy_map[query_kmer_exact]
                                            [target_kmer_fuzzy]):
                        query_kmer_alignments['+'].append(
                            [target_kmer_pos, query_kmer_pos, False])
        # query reverse complement alignment
        if args.seq_type == 'dna':
            query_kmer_rc = reverse_complement(query_kmer)
            query_kmer_rc_exact = ''.join([query_kmer_rc[x]
                                           for x in args.seed_exact_idxs])
            query_kmer_rc_fuzzy = ''.join([query_kmer_rc[x]
                                           for x in args.seed_fuzzy_idxs])
            if query_kmer_rc_exact in fuzzy_map:
                for target_kmer_fuzzy in fuzzy_map[query_kmer_rc_exact]:
                    similarity_score = calc_similarity_score(
                        query_kmer_rc_fuzzy, target_kmer_fuzzy,
                        args.seq_type, args.sim_algo, aligner)
                    if similarity_score >= args.sim_cutoff:
                        for target_kmer_pos in (fuzzy_map[query_kmer_rc_exact]
                                                [target_kmer_fuzzy]):
                            query_kmer_alignments['-'].append(
                                [target_kmer_pos, query_kmer_pos, False])
    # group query k-mer alignments
    cdef dict alignment
    cdef str strand
    cdef list alignments, alignment_grp, alignment_grp_idxs
    cdef np.ndarray alignment_grp_scores = np.zeros(
        (len(query_seq),), dtype=np.bool)
    cdef unsigned int i, j, grp_target_seq_pos, grp_query_seq_pos
    for strand, alignments in query_kmer_alignments.items():
        alignments.sort(key=itemgetter(0, 1))
        for i in range(len(alignments)):
            if not alignments[i][2]:
                alignment_grp = [tuple(alignments[i][:2])]
                alignment_grp_idxs = [i]
                for j in range(i + 1, len(alignments)):
                    if (alignments[j][0] - (alignment_grp[-1][0] + args.k) <=
                            args.max_kmer_gap):
                        if ((alignments[j][0] < alignment_grp[-1][0] + args.k)
                                or (alignments[j][1] < alignment_grp[-1][1]
                                    + args.k)):
                            if (alignments[j][0] - alignment_grp[-1][0] ==
                                    alignments[j][1] - alignment_grp[-1][1]):
                                alignment_grp.append(tuple(alignments[j][:2]))
                                alignment_grp_idxs.append(j)
                        elif (0 < alignments[j][1] - (alignment_grp[-1][1]
                                                      + args.k) <=
                              args.max_kmer_gap):
                            alignment_grp.append(tuple(alignments[j][:2]))
                            alignment_grp_idxs.append(j)
                    else:
                        # flag alignments used in group
                        for idx in alignment_grp_idxs:
                            alignments[idx][2] = True
                        alignment = build_alignment(alignment_grp, query_seq,
                                                    target_seq, aligner, args)
                        if alignment['e_value'] <= args.expect_thres:
                            yield alignment
                        break
