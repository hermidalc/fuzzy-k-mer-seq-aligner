# Fuzzy K-mer Local Sequence Aligner

```bash
$ ./fuzzy_kmer_seq_aligner.py --help
usage: fuzzy_kmer_seq_aligner.py [-h] -fs FUZZY_SEED
                                 -st {dna,rna,protein}
                                 -qs QUERY_SEQ
                                 -ts TARGET_SEQ
                                 [-ql QUERY_LOC QUERY_LOC]
                                 [-tl TARGET_LOC TARGET_LOC]
                                 [-sa {levenshtein,hamming,smith-waterman}]
                                 [-sc SIM_CUTOFF]
                                 [-ms MATCH_SCORE]
                                 [-ss MISMATCH_SCORE]
                                 [-og OPEN_GAP_SCORE]
                                 [-eg EXTEND_GAP_SCORE]
                                 [-sm SUB_MATRIX]
                                 [-mg MAX_KMER_GAP]
                                 [-et EXPECT_THRES]
                                 [-af {pairwise,tabular}]
                                 [-as {e_value,bit_score,pct_id,q_start,s_start}]
                                 [-ma MAX_ALIGNS]
                                 [-pw PWA_WIDTH]

optional arguments:
  -h, --help            show this help message and exit
  -fs FUZZY_SEED, --fuzzy-seed FUZZY_SEED
                        fuzzy k-mer seed pattern
  -st {dna,rna,protein}, --seq-type {dna,rna,protein}
                        sequence type
  -qs QUERY_SEQ, --query-seq QUERY_SEQ
                        query FASTA sequence file
  -ts TARGET_SEQ, --target-seq TARGET_SEQ
                        target FASTA sequence file
  -ql QUERY_LOC QUERY_LOC, --query-loc QUERY_LOC QUERY_LOC
                        query sequence start stop
  -tl TARGET_LOC TARGET_LOC, --target-loc TARGET_LOC TARGET_LOC
                        target sequence start stop
  -sa {levenshtein,hamming,smith-waterman}, --sim-algo {levenshtein,hamming,smith-waterman}
                        string similarity algorithm
  -sc SIM_CUTOFF, --sim-cutoff SIM_CUTOFF
                        fuzzy membership similarity cutoff
  -ms MATCH_SCORE, --match-score MATCH_SCORE
                        match score
  -ss MISMATCH_SCORE, --mismatch-score MISMATCH_SCORE
                        mismatch score
  -og OPEN_GAP_SCORE, --open-gap-score OPEN_GAP_SCORE
                        open gap score
  -eg EXTEND_GAP_SCORE, --extend-gap-score EXTEND_GAP_SCORE
                        extend gap score
  -sm SUB_MATRIX, --sub-matrix SUB_MATRIX
                        substitution matrix (any Biopython MatrixInfo matrix
                        name)
  -mg MAX_KMER_GAP, --max-kmer-gap MAX_KMER_GAP
                        maximum gap size when grouping kmers
  -et EXPECT_THRES, --expect-thres EXPECT_THRES
                        expect value threshold
  -af {pairwise,tabular}, --align-fmt {pairwise,tabular}
                        alignment output format
  -as {e_value,bit_score,pct_id,q_start,s_start}, --align-sort {e_value,bit_score,pct_id,q_start,s_start}
                        alignment output sort
  -ma MAX_ALIGNS, --max-aligns MAX_ALIGNS
                        maximum number of alignments to show
  -pw PWA_WIDTH, --pwa-width PWA_WIDTH
                        pairwise alignment output width
```

## Installation

Install latest
<a href="https://www.anaconda.com/distribution/" target="_blank">
Anaconda
</a>
or
<a href="https://docs.conda.io/en/latest/miniconda.html" target="_blank">
Miniconda
</a>
for Python 3.x.

Update ``conda`` to latest:
```bash
conda update -n base conda
```

Clone the Github repository:
```bash
git clone git@github.com:hermidalc/fuzzy-kmer-seq-aligner.git
```

Create and activate the conda environment:
```bash
cd fuzzy-kmer-seq-aligner
conda env create -f environment.yml
conda activate fuzzy-kmer-seq-aligner
```

Part of the fuzzy k-mer aligner program is written in Cython and this
extension gets automatically compiled when you first run it. Under
Linux no additional dependencies need to be installed. Under Mac
OSX you will need the Apple Clang Compiler tools for Cython to work.
Download and install
<a href="https://developer.apple.com/download/">
Apple Xcode
</a>
to get the compiler tools. Under Windows you will need the Microsoft
Visual C++ compiler tools.  The version must be the same as the MSVC++
compiler used to build the Python executable you are using. For Python
3.7 download and install
<a href="https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2019" target="_blank">
Build Tools for Visual Studio 2019
</a>.

In order to run all ``snakemake`` examples the BLAST+ suite must also be
installed. Under Linux and Mac OSX install ``blast`` in your activated
conda environment:
```bash
conda install blast
```

Conda has no BLAST+ suite package for Windows.  Download and install the
BLAST+ suite from the NCBI FTP site (see
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST).  The installer should
automatically add the BLAST executable directories to your ``%PATH%``.

## Usage

To run all examples simply type ``snakemake all`` and to clean up example
data and results type ``snakemake clean``.
