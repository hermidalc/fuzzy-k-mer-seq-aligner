# functions.pyx
# cython: language_level=3, boundscheck=False, wraparound=True, cdivision=True

import numpy as np
cimport numpy as np
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
def build_fuzzy_map(str seq, unsigned int k, list seed_exact_idxs,
                    list seed_fuzzy_idxs):
    cdef dict fuzzy_map = {}
    cdef unsigned int kmer_pos
    cdef str kmer, kmer_exact, kmer_fuzzy
    for kmer_pos, kmer in kmers(seq, k):
        kmer_exact = ''.join([kmer[x] for x in seed_exact_idxs])
        kmer_fuzzy = ''.join([kmer[x] for x in seed_fuzzy_idxs])
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
        score = 0.
    return score


# build alignments from query k-mer alignments
def build_alignments(dict fuzzy_map, str query_seq, str target_seq,
                     unsigned int k, list seed_exact_idxs,
                     list seed_fuzzy_idxs, aligner, args):
    # align query k-mers to target sequence
    cdef dict query_kmer_alignments = {'+': [], '-': []}
    cdef unsigned int query_kmer_pos, target_kmer_pos
    cdef float similarity_score
    cdef str query_kmer, query_kmer_rc, query_kmer_exact, query_kmer_rc_exact
    cdef str query_kmer_fuzzy, query_kmer_rc_fuzzy, target_kmer_fuzzy
    cdef np.ndarray query_kmers = np.empty(
        (len(query_seq) - k + 1,), dtype=np.dtype('U' + str(k)))
    for query_kmer_pos, query_kmer in kmers(query_seq, k):
        query_kmers[query_kmer_pos] = query_kmer
        query_kmer_exact = ''.join([query_kmer[x] for x in seed_exact_idxs])
        query_kmer_fuzzy = ''.join([query_kmer[x] for x in seed_fuzzy_idxs])
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
                                           for x in seed_exact_idxs])
            query_kmer_rc_fuzzy = ''.join([query_kmer_rc[x]
                                           for x in seed_fuzzy_idxs])
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
    # build query alignments from k-mer alignments
    cdef str strand
    cdef list alignments, alignment_grp
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
                    if alignments[j][0] - alignment_grp[-1][0] < k:
                        if 0 < alignments[j][1] - alignment_grp[-1][1] < k:
                            alignment_grp.append(tuple(alignments[j][:2]))
                            alignment_grp_idxs.append(j)
                    else:
                        for idx in alignment_grp_idxs:
                            alignments[idx][2] = True
                        # TODO: build alignments from alignment group here
                        break
