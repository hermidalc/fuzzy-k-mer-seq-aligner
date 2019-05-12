# Fuzzy K-mer Local Sequence Aligner

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
Under Windows you must install BLAST+ manually from the
<a href="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST" target="_blank">NCBI BLAST FTP</a> site.
