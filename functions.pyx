# functions.pyx
# cython: language_level=3, boundscheck=False, wraparound=True, cdivision=True

from math import log
import numpy as np
cimport numpy as np
from operator import itemgetter
from Bio.Seq import reverse_complement
from Levenshtein import distance, hamming
from ka_config import KA_NA_ROUND_DOWN


# iterate over kmers in a sequence
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
                          dict aligners):
    cdef str c
    cdef float score = 0
    cdef float perfect_score
    if algo == 'levenshtein':
        score = 1 - distance(query, target) / len(target)
    elif algo == 'hamming':
        score = 1 - hamming(query, target) / len(target)
    elif algo == 'smith-waterman':
        if seq_type in ('dna', 'rna'):
            perfect_score = len(target) * aligners['local'].match
        else:
            perfect_score = np.sum(
                [aligners['local'].substitution_matrix[(c, c)]
                for c in target])
        score = aligners['local'].score(query, target) / perfect_score
    return score


# build local alignment from query kmer alignment group
def build_alignment(list alignment_grp, str query_seq, str target_seq,
                    str strand, dict aligners, args):
    cdef list query_align_parts = []
    cdef list match_align_parts = []
    cdef list target_align_parts = []
    cdef list match_align_gap_chrs
    cdef str query_align_gap, match_align_gap, target_align_gap
    cdef str query_gap, query_kmer, query_chr
    cdef str target_gap, target_kmer, target_chr
    cdef str query_alignment, match_alignment, target_alignment
    cdef str c
    cdef unsigned int query_seq_pos, query_seq_start_pos
    cdef unsigned int target_seq_pos, target_seq_start_pos
    cdef unsigned int num_idents = 0
    cdef unsigned int num_gaps = 0
    cdef unsigned int i
    cdef int query_gap_size, target_gap_size, g
    cdef int raw_score = 0
    cdef int query_seq_curr_pos = -1
    cdef int target_seq_curr_pos = -1
    cdef float bit_score
    for target_seq_pos, query_seq_pos in alignment_grp:
        if target_seq_curr_pos >= 0:
            target_gap_size = target_seq_pos - target_seq_curr_pos
            query_gap_size = query_seq_pos - query_seq_curr_pos
            if target_gap_size < 0:
                target_seq_start_pos = target_seq_curr_pos
                query_seq_start_pos = query_seq_curr_pos
            elif target_gap_size > 0:
                if strand == 'Plus':
                    target_gap = target_seq[
                        target_seq_curr_pos:target_seq_pos]
                else:
                    target_gap = reverse_complement(target_seq[
                        (len(target_seq) - target_seq_pos):
                        (len(target_seq) - target_seq_curr_pos)])
                if query_gap_size > 0:
                    query_gap = query_seq[query_seq_curr_pos:(
                        query_seq_curr_pos + query_gap_size)]
                    gap_alignment = aligners['global'].align(
                        query_gap, target_gap)[0]
                    raw_score += gap_alignment.score
                    query_align_gap, match_align_gap, target_align_gap = str(
                        gap_alignment).split('\n')[:3]
                    if args.seq_type == 'protein':
                        match_align_gap_chrs = list(match_align_gap)
                        for i, c in enumerate(match_align_gap_chrs):
                            if c == '|':
                                match_align_gap_chrs[i] = query_align_gap[i]
                            elif c == 'X':
                                if (aligners['global'].substitution_matrix[
                                        (query_align_gap[i],
                                         target_align_gap[i]) ] > 0):
                                    match_align_gap_chrs[i] = '+'
                                else:
                                    match_align_gap_chrs[i] = ' '
                            elif c == '-':
                                match_align_gap_chrs[i] = ' '
                        match_align_gap = ''.join(match_align_gap_chrs)
                    query_align_parts.append(query_align_gap)
                    match_align_parts.append(match_align_gap)
                    target_align_parts.append(target_align_gap)
                else:
                    raw_score += sum(
                        args.open_gap_score if g == 0 else
                        args.extend_gap_score for g in range(target_gap_size))
                    query_align_parts.append('-' * target_gap_size)
                    match_align_parts.append(' ' * target_gap_size)
                    target_align_parts.append(target_gap)
                target_seq_start_pos = target_seq_pos
                query_seq_start_pos = query_seq_pos
            else:
                if query_gap_size > 0:
                    query_gap = query_seq[query_seq_curr_pos:(
                        query_seq_curr_pos + query_gap_size)]
                    raw_score += sum(
                        args.open_gap_score if g == 0 else
                        args.extend_gap_score for g in range(query_gap_size))
                    query_align_parts.append(query_gap)
                    match_align_parts.append(' ' * query_gap_size)
                    target_align_parts.append('-' * query_gap_size)
                target_seq_start_pos = target_seq_pos
                query_seq_start_pos = query_seq_pos
        else:
            target_seq_start_pos = target_seq_pos
            query_seq_start_pos = query_seq_pos
        target_seq_curr_pos = target_seq_pos + args.k
        query_seq_curr_pos = query_seq_pos + args.k
        if strand == 'Plus':
            target_kmer = target_seq[
                target_seq_start_pos:target_seq_curr_pos]
        else:
            target_kmer = reverse_complement(target_seq[
                (len(target_seq) - target_seq_curr_pos):
                (len(target_seq) - target_seq_start_pos)])
        query_kmer = query_seq[query_seq_start_pos:query_seq_curr_pos]
        query_align_parts.append(query_kmer)
        target_align_parts.append(target_kmer)
        for target_chr, query_chr in zip(target_kmer, query_kmer):
            if args.seq_type in ('dna', 'rna'):
                if query_chr == target_chr:
                    raw_score += args.match_score
                    match_align_parts.append('|')
                else:
                    raw_score += args.mismatch_score
                    match_align_parts.append(' ')
            else:
                raw_score += aligners['global'].substitution_matrix[
                    (query_chr, target_chr)]
                if query_chr == target_chr:
                    match_align_parts.append(query_chr)
                elif (aligners['global'].substitution_matrix[
                        (query_chr, target_chr)] > 0):
                    match_align_parts.append('+')
                else:
                    match_align_parts.append(' ')
    query_alignment = ''.join(query_align_parts)
    match_alignment = ''.join(match_align_parts)
    target_alignment = ''.join(target_align_parts)
    if strand == 'Plus':
        tstart = alignment_grp[0][0] + 1
        tend = alignment_grp[-1][0] + args.k
    else:
        tstart = len(target_seq) - alignment_grp[0][0]
        tend = len(target_seq) - (alignment_grp[-1][0] + args.k) + 1
    if (args.seq_type in ('dna', 'rna') and raw_score % 2 != 0 and
            (args.match_score, args.mismatch_score) in KA_NA_ROUND_DOWN):
        raw_score = max(raw_score - 1, 0)
    bit_score = (args.ka_gapped_l * raw_score - log(args.ka_gapped_k)) / log(2)
    if args.seq_type in ('dna', 'rna'):
        num_ids = match_alignment.count('|')
    else:
        num_ids = np.sum([match_alignment.count(c)
                          for c in args.seq_abc.letters])
    return {'query': query_alignment,
            'match': match_alignment,
            'target': target_alignment,
            'qstart': alignment_grp[0][1] + 1,
            'qend': alignment_grp[-1][1] + args.k,
            'tstart': tstart,
            'tend': tend,
            'raw_score': raw_score,
            'bit_score': bit_score,
            'e_value': len(query_seq) * len(target_seq) / (2 ** bit_score),
            'p_value': 2 ** (-bit_score),
            'num_ids': num_ids,
            'num_pos': match_alignment.count('+'),
            'num_gaps': (query_alignment.count('-') +
                         target_alignment.count('-')),
            'strand': strand}


def pairwise_align(dict fuzzy_map, str query_seq, str target_seq,
                   dict aligners, args):
    # align query kmers to target sequence
    cdef dict query_kmer_alignments = {'Plus': [], 'Minus': []}
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
                if target_kmer_fuzzy != '':
                    similarity_score = calc_similarity_score(
                        query_kmer_fuzzy, target_kmer_fuzzy, args.seq_type,
                        args.sim_algo, aligners)
                else:
                    similarity_score = 1
                if similarity_score >= args.sim_cutoff:
                    for target_kmer_pos in (fuzzy_map[query_kmer_exact]
                                            [target_kmer_fuzzy]):
                        query_kmer_alignments['Plus'].append(
                            [target_kmer_pos, query_kmer_pos, False])
        # align query kmer reverse complement
        if args.seq_type == 'dna':
            query_kmer_rc = reverse_complement(query_kmer)
            query_kmer_rc_exact = ''.join([query_kmer_rc[x]
                                           for x in args.seed_exact_idxs])
            query_kmer_rc_fuzzy = ''.join([query_kmer_rc[x]
                                           for x in args.seed_fuzzy_idxs])
            if query_kmer_rc_exact in fuzzy_map:
                for target_kmer_fuzzy in fuzzy_map[query_kmer_rc_exact]:
                    if target_kmer_fuzzy != '':
                        similarity_score = calc_similarity_score(
                            query_kmer_rc_fuzzy, target_kmer_fuzzy,
                            args.seq_type, args.sim_algo, aligners)
                    else:
                        similarity_score = 1
                    if similarity_score >= args.sim_cutoff:
                        for target_kmer_pos in (fuzzy_map[query_kmer_rc_exact]
                                                [target_kmer_fuzzy]):
                            query_kmer_alignments['Minus'].append(
                                [len(target_seq) - (target_kmer_pos + args.k),
                                 query_kmer_pos, False])
    # group query kmer alignments
    cdef dict alignment
    cdef str strand
    cdef list alignments, alignment_grp, alignment_grp_idxs
    cdef unsigned int i, j, k
    for strand, alignments in query_kmer_alignments.items():
        alignments.sort(key=itemgetter(0, 1))
        for i in range(len(alignments)):
            if not alignments[i][2]:
                alignment_grp = [tuple(alignments[i][:2])]
                alignment_grp_idxs = [i]
                for j in range(i + 1, len(alignments)):
                    if (alignments[j][0] - (alignment_grp[-1][0] + args.k)
                            <= args.max_kmer_gap):
                        if ((alignments[j][0] < alignment_grp[-1][0] + args.k)
                                or (alignments[j][1] < alignment_grp[-1][1]
                                    + args.k)):
                            if (alignments[j][0] - alignment_grp[-1][0] ==
                                    alignments[j][1] - alignment_grp[-1][1]):
                                alignment_grp.append(tuple(alignments[j][:2]))
                                alignment_grp_idxs.append(j)
                        elif (0 < alignments[j][1] -
                              (alignment_grp[-1][1] + args.k)
                              <= args.max_kmer_gap):
                            alignment_grp.append(tuple(alignments[j][:2]))
                            alignment_grp_idxs.append(j)
                    else:
                        break
                alignment = build_alignment(alignment_grp, query_seq,
                                            target_seq, strand, aligners,
                                            args)
                if alignment['e_value'] <= args.expect_thres:
                    yield alignment
                # flag alignments in group as used
                for k in alignment_grp_idxs:
                    alignments[k][2] = True
