# Fuzzy K-mer Local Sequence Aligner

```bash
$ ./fuzzy_kmer_seq_aligner.py --help
usage: fuzzy_kmer_seq_aligner.py [-h] --fuzzy-seed FUZZY_SEED
                                 --query-seq QUERY_SEQ
                                 [--query-loc QUERY_LOC QUERY_LOC]
                                 --target-seq TARGET_SEQ
                                 [--target-loc TARGET_LOC TARGET_LOC]
                                 --seq-type {dna,rna,protein}
                                 [--sim-algo {levenshtein,hamming,smith-waterman}]
                                 [--sim-cutoff SIM_CUTOFF]
                                 [--match-score MATCH_SCORE]
                                 [--mismatch-score MISMATCH_SCORE]
                                 [--open-gap-score OPEN_GAP_SCORE]
                                 [--extend-gap-score EXTEND_GAP_SCORE]
                                 [--sub-matrix SUB_MATRIX]
                                 [--max-kmer-gap MAX_KMER_GAP]
                                 [--expect-thres EXPECT_THRES]
                                 [--align-fmt {pairwise,tabular}]
                                 [--align-sort {e_value,bit_score,pct_id,q_start,s_start}]
                                 [--max-aligns MAX_ALIGNS]
                                 [--pwa-width PWA_WIDTH]

optional arguments:
  -h, --help            show this help message and exit
  --fuzzy-seed FUZZY_SEED, -fs FUZZY_SEED
                        fuzzy k-mer seed pattern
  --query-seq QUERY_SEQ, -qs QUERY_SEQ
                        query FASTA sequence file
  --query-loc QUERY_LOC QUERY_LOC, -ql QUERY_LOC QUERY_LOC
                        query sequence start-stop
  --target-seq TARGET_SEQ, -ts TARGET_SEQ
                        target FASTA sequence file
  --target-loc TARGET_LOC TARGET_LOC, -tl TARGET_LOC TARGET_LOC
                        target sequence start-stop
  --seq-type {dna,rna,protein}, -st {dna,rna,protein}
                        sequence type
  --sim-algo {levenshtein,hamming,smith-waterman}, -sa {levenshtein,hamming,smith-waterman}
                        string similarity algorithm
  --sim-cutoff SIM_CUTOFF, -sc SIM_CUTOFF
                        fuzzy membership similarity cutoff
  --match-score MATCH_SCORE, -ms MATCH_SCORE
                        match score
  --mismatch-score MISMATCH_SCORE, -ss MISMATCH_SCORE
                        mismatch score
  --open-gap-score OPEN_GAP_SCORE, -og OPEN_GAP_SCORE
                        open gap score
  --extend-gap-score EXTEND_GAP_SCORE, -eg EXTEND_GAP_SCORE
                        extend gap score
  --sub-matrix SUB_MATRIX, -sm SUB_MATRIX
                        substitution matrix (any Biopython MatrixInfo matrix
                        name)
  --max-kmer-gap MAX_KMER_GAP, -mg MAX_KMER_GAP
                        maximum gap size when grouping kmers
  --expect-thres EXPECT_THRES, -et EXPECT_THRES
                        expect value threshold
  --align-fmt {tabular,pairwise}, -af {tabular,pairwise}
                        alignment output format
  --align-sort {e_value,bit_score,pct_id,q_start,s_start}, -as {e_value,bit_score,pct_id,q_start,s_start}
                        alignment output sort
  --max-aligns MAX_ALIGNS, -ma MAX_ALIGNS
                        maximum number of alignments to show
  --pwa-width PWA_WIDTH, -aw PWA_WIDTH
                        pairwise alignment output width
```

## Installation

Install latest
<a href="https://www.anaconda.com/distribution/" target="_blank">Anaconda</a>
or
<a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">Miniconda</a>
for Python 3.x.

Update ``conda`` to latest:
```bash
conda update -n base conda
```

Clone the Github repository:
```bash
git clone git@github.com:hermidalc/fuzzy-kmer-seq-aligner.git
```

Create and activate conda environment:
```bash
conda env create -f environment.yml
conda activate fuzzy-kmer-seq-aligner
```

## Running examples

In order to run all ``snakemake`` examples the BLAST+ suite must be installed.
Under Linux and Mac OSX install ``blast`` in your activated conda environment:
```bash
conda install blast
```
Under Windows you must manually download and install the BLAST+ suite from the
NCBI FTP site (see ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST) and
ensure that the BLAST executables are in your ``%PATH%``.

To run all examples simply type ``snakemake all`` and to clean up example
data and results type ``snakemake clean``
