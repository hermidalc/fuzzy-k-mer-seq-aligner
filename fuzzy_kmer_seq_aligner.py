#!/usr/bin/env python

from argparse import ArgumentParser
from importlib import import_module
from os import path
from textwrap import dedent
from Bio.Align import PairwiseAligner
from Bio.Alphabet.IUPAC import (ExtendedIUPACProtein, IUPACAmbiguousDNA,
                                IUPACAmbiguousRNA)
from Bio.SeqIO.FastaIO import SimpleFastaParser
from ka_config import KA_PARAMS
import pyximport
pyximport.install(inplace=True)
from functions import build_fuzzy_map, pairwise_align

parser = ArgumentParser()
parser.add_argument('--fuzzy-seed', '-fs', type=str, required=True,
                    help='fuzzy k-mer seed pattern')
parser.add_argument('--query-seq', '-qs', type=str, required=True,
                    help='query FASTA sequence file')
parser.add_argument('--target-seq', '-ts', type=str, required=True,
                    help='target FASTA sequence file')
parser.add_argument('--seq-type', '-st', type=str, required=True,
                    choices=['dna', 'rna', 'protein'], help='sequence type')
parser.add_argument('--sim-algo', '-sa', type=str, default='smith-waterman',
                    choices=['levenshtein', 'hamming', 'smith-waterman'],
                    help='string similarity algorithm')
parser.add_argument('--sim-cutoff', '-sc', type=float, default='0.6',
                    help='fuzzy membership similarity cutoff')
parser.add_argument('--match-score', '-ms', type=float, default='2',
                    help='match score')
parser.add_argument('--mismatch-score', '-ss', type=float, default='-3',
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
parser.add_argument('--expect-thres', '-et', type=float, default='10',
                    help='expect value threshold')
parser.add_argument('--align-fmt', '-af', type=str, default='pairwise',
                    choices=['tabular', 'pairwise'],
                    help='alignment output format')
parser.add_argument('--align-sort', '-as', type=str, default='e_value',
                    choices=['e_value', 'score', 'pct_id', 'q_start',
                             's_start'],
                    help='alignment output sort')
parser.add_argument('--max-aligns', '-ma', type=int, default='50',
                    help='maximum number of alignments to show')
parser.add_argument('--pwa-width', '-aw', type=int, default='60',
                    help='pairwise alignment output width')
args = parser.parse_args()

for file in (args.query_seq, args.target_seq):
    if not path.isfile(file):
        parser.error("File %s doesn't exist" % file)

# alphabet
if args.seq_type == 'dna':
    seq_alphabet = IUPACAmbiguousDNA()
elif args.seq_type == 'rna':
    seq_alphabet = IUPACAmbiguousRNA()
else:
    seq_alphabet = ExtendedIUPACProtein()

# Aligners setup
aligners = {'global': PairwiseAligner(), 'local': None}
aligners['global'].mode = 'global'
if args.seq_type in ('dna', 'rna'):
    aligners['global'].match = args.match_score
    aligners['global'].mismatch = args.mismatch_score
    if not args.open_gap_score:
        args.open_gap_score = -5
    if not args.extend_gap_score:
        args.extend_gap_score = -2
else:
    sub_matrix = getattr(import_module('Bio.SubsMat.MatrixInfo'),
                         args.sub_matrix)
    aligners['global'].substitution_matrix = sub_matrix
    if not args.open_gap_score:
        args.open_gap_score = -11
    if not args.extend_gap_score:
        args.extend_gap_score = -1
aligners['global'].open_gap_score = args.open_gap_score
aligners['global'].extend_gap_score = args.extend_gap_score
if args.sim_algo == 'smith-waterman':
    aligners['local'] = PairwiseAligner()
    aligners['local'].mode = 'local'
    if args.seq_type in ('dna', 'rna'):
        aligners['local'].match = args.match_score
        aligners['local'].mismatch = args.mismatch_score
    else:
        aligners['local'].substitution_matrix = sub_matrix
    aligners['local'].open_gap_score = args.open_gap_score
    aligners['local'].extend_gap_score = args.extend_gap_score

# Karlin-Altschul parameter values
if args.seq_type in ('dna', 'rna'):
    if ((args.match_score, args.mismatch_score) in KA_PARAMS['na'] and
            (abs(args.open_gap_score), abs(args.extend_gap_score)) in
            KA_PARAMS['na'][(args.match_score, args.mismatch_score)]):
        args.ka_gapped_l = KA_PARAMS['na'][
            (args.match_score, args.mismatch_score)][
                (abs(args.open_gap_score), abs(args.extend_gap_score))][0]
        args.ka_gapped_k = KA_PARAMS['na'][
            (args.match_score, args.mismatch_score)][
                (abs(args.open_gap_score), abs(args.extend_gap_score))][1]
    else:
        args.ka_gapped_l = 1.280
        args.ka_gapped_k = 0.460
else:
    if (args.sub_matrix in KA_PARAMS['aa'] and
            (abs(args.open_gap_score), abs(args.extend_gap_score)) in
            KA_PARAMS['aa'][args.sub_matrix]):
        args.ka_gapped_l = KA_PARAMS['aa'][args.sub_matrix][
            (abs(args.open_gap_score), abs(args.extend_gap_score))][0]
        args.ka_gapped_k = KA_PARAMS['aa'][args.sub_matrix][
            (abs(args.open_gap_score), abs(args.extend_gap_score))][1]
    else:
        args.ka_gapped_l = 0.267
        args.ka_gapped_k = 0.041

if args.align_sort in ('e_value', 'q_start', 's_start'):
    align_sort_rev = False
elif args.align_sort in ('score', 'pct_id'):
    align_sort_rev = True

if args.align_fmt == 'pairwise':
    pw_header_fmt = dedent('''\
    Query: {qid} (Length = {qlen})
    Sbjct: {sid} (Length = {slen})
    ''')
    pw_section_header_fmt = dedent('''\
    Score = {bits:.1f} bits ({raw}), Expect = {eval:{efmt}}
    Identities = {ids}/{idt} ({idp:.0%}), Gaps = {gaps}/{gapt} ({gapp:.0%})\
    {strand}
    ''')
    pw_alignment_fmt = dedent('''\
    Query   {qstar:<{lenpad}}   {query}   {qend}
    {mpad}{match}
    Sbjct   {sstar:<{lenpad}}   {sbjct}   {send}
    ''')

args.k = len(args.fuzzy_seed)
args.seed_exact_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '#']
args.seed_fuzzy_idxs = [i for i, c in enumerate(args.fuzzy_seed) if c == '*']
query_seq_fh = open(args.query_seq, 'r')
target_seq_fh = open(args.target_seq, 'r')
for target_seq_title, target_seq in SimpleFastaParser(target_seq_fh):
    fuzzy_map = build_fuzzy_map(target_seq, args)
    for query_seq_title, query_seq in SimpleFastaParser(query_seq_fh):
        if args.align_fmt == 'pairwise':
            lenpad = len(str(max(len(query_seq), len(target_seq))))
            print(pw_header_fmt.format(
                qid=query_seq_title, qlen=len(query_seq),
                sid=target_seq_title, slen=len(target_seq)))
        for i, alignment in enumerate(sorted(
                pairwise_align(fuzzy_map, query_seq, target_seq, aligners,
                               args),
                key=lambda a: a[args.align_sort], reverse=align_sort_rev)):
            if args.align_fmt == 'pairwise':
                efmt = '.3' + ('e' if alignment['e_value'] < 1e-3 else 'f')
                idp = alignment['num_ids'] / len(query_seq) * 100
                gapp = alignment['num_gaps'] / len(query_seq) * 100
                if args.seq_type == 'dna':
                    strand = '\nStrand = Plus/' + alignment['strand']
                else:
                    strand = ''
                print(pw_section_header_fmt.format(
                    bits=alignment['bit_score'], raw=alignment['raw_score'],
                    eval=alignment['e_value'], efmt=efmt,
                    ids=alignment['num_ids'], idt=len(query_seq), idp=idp,
                    gaps=alignment['num_gaps'], gapt=len(query_seq),
                    gapp=gapp, strand=strand))
                mpad = ' ' * (8 + lenpad + 3)
                for q, m, t in zip(
                        range(alignment['qstart'], alignment['qend'] + 1,
                              args.pwa_width),
                        range(0, len(alignment['match']), args.pwa_width),
                        range(alignment['tstart'], alignment['tend'] + 1,
                              args.pwa_width)):
                    if m + args.pwa_width < len(alignment['match']):
                        pwa_width = args.pwa_width
                    else:
                        pwa_width = len(alignment['match'])
                    if alignment['strand'] == 'Plus':
                        sstar = t + 1
                        send = t + pwa_width
                    else:
                        sstar = t + pwa_width
                        send = t + 1
                    print(pw_alignment_fmt.format(
                        qstar=q + 1, lenpad=lenpad,
                        query=alignment['query'][m:(m + pwa_width)],
                        qend=q + pwa_width,
                        match=alignment['match'][m:(m + pwa_width)],
                        sstar=sstar, mpad=mpad,
                        sbjct=alignment['target'][m:(m + pwa_width)],
                        send=send))
            if i + 1 == args.max_aligns:
                break
    del fuzzy_map
target_seq_fh.close()
query_seq_fh.close()
