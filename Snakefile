################################################################################
# SETUP
################################################################################
# Modules
from os.path import join
import shutil

# Configuration
OUTPUT_DIR = config['output_dir'] = config.get('output_dir', 'results')

# Directories
DATA_DIR = 'data'
SRC_DIR = '.'

# Files
EX1_QUERY_DATA = join(DATA_DIR, 'NG_007458.fna')
EX1_TARGET_DATA = join(DATA_DIR, 'NC_000068.fna')
EX2_QUERY_DATA = join(DATA_DIR, 'Q9S3R8.faa')
EX2_TARGET_DATA = join(DATA_DIR, 'Q9S3R9.faa')

EX1_FUZZY_PW_ALIGN_OUT = join(OUTPUT_DIR, 'ex1_fuzzy_pw_align.txt')
EX1_BLAST_PW_ALIGN_OUT = join(OUTPUT_DIR, 'ex1_blast_pw_align.txt')
EX2_FUZZY_PW_ALIGN_OUT = join(OUTPUT_DIR, 'ex2_fuzzy_pw_align.txt')
EX2_BLAST_PW_ALIGN_OUT = join(OUTPUT_DIR, 'ex2_blast_pw_align.txt')

# Scripts
FUZZY_KMER_ALIGNER = join(SRC_DIR, 'fuzzy_kmer_seq_aligner.py')
BLASTN = 'blastn'
BLASTP = 'blastp'

################################################################################
# GENRAL RULES
################################################################################

rule all:
    input:
        EX1_QUERY_DATA,
        EX1_TARGET_DATA,
        EX1_FUZZY_PW_ALIGN_OUT,
        EX1_BLAST_PW_ALIGN_OUT,
        EX2_QUERY_DATA,
        EX2_TARGET_DATA,
        EX2_FUZZY_PW_ALIGN_OUT,
        EX2_BLAST_PW_ALIGN_OUT

rule clean:
    run:
        shutil.rmtree(DATA_DIR, ignore_errors=True)
        shutil.rmtree(OUTPUT_DIR, ignore_errors=True)

rule get_ex1_query_data:
    params:
        url="https://www.ncbi.nlm.nih.gov/search/api/sequence/NG_007458/?report=fasta"
    output:
        EX1_QUERY_DATA
    shell:
        "wget -O {output} '{params.url}'"

rule get_ex1_target_data:
    params:
        url="https://www.ncbi.nlm.nih.gov/search/api/sequence/NC_000068/?report=fasta&from=26400000&to=26599999"
    output:
        EX1_TARGET_DATA
    shell:
        "wget -O {output} '{params.url}'"

rule run_ex1_fuzzy_kmers:
    input:
        EX1_QUERY_DATA,
        EX1_TARGET_DATA
    output:
        EX1_FUZZY_PW_ALIGN_OUT
    shell:
        "{FUZZY_KMER_ALIGNER} -st dna -qs {EX1_QUERY_DATA} "
        "-ts {EX1_TARGET_DATA} -fs '########*************' "
        "-as bit_score > {EX1_FUZZY_PW_ALIGN_OUT}"

rule run_ex1_blast:
    input:
        EX1_QUERY_DATA,
        EX1_TARGET_DATA
    output:
        EX1_BLAST_PW_ALIGN_OUT
    shell:
        "{BLASTN} -task blastn -query {EX1_QUERY_DATA} "
        "-subject {EX1_TARGET_DATA} -out {EX1_BLAST_PW_ALIGN_OUT} "
        "-word_size 11 -reward 2 -penalty -3 -gapopen 5 -gapextend 2 "
        "-outfmt 0 -dust no -soft_masking false"

rule get_ex2_query_data:
    params:
        url="https://www.ncbi.nlm.nih.gov/search/api/sequence/Q9S3R8/?report=fasta"
    output:
        EX2_QUERY_DATA
    shell:
        "wget -O {output} '{params.url}'"

rule get_ex2_target_data:
    params:
        url="https://www.ncbi.nlm.nih.gov/search/api/sequence/Q9S3R9/?report=fasta"
    output:
        EX2_TARGET_DATA
    shell:
        "wget -O {output} '{params.url}'"

rule run_ex2_fuzzy_kmers:
    input:
        EX2_QUERY_DATA,
        EX2_TARGET_DATA
    output:
        EX2_FUZZY_PW_ALIGN_OUT
    shell:
        "{FUZZY_KMER_ALIGNER} -st protein -qs {EX2_QUERY_DATA} "
        "-ts {EX2_TARGET_DATA} -fs '###*****' "
        "-as bit_score > {EX2_FUZZY_PW_ALIGN_OUT}"

rule run_ex2_blast:
    input:
        EX2_QUERY_DATA,
        EX2_TARGET_DATA
    output:
        EX2_BLAST_PW_ALIGN_OUT
    shell:
        "{BLASTP} -task blastp -query {EX2_QUERY_DATA} "
        "-subject {EX2_TARGET_DATA} -out {EX2_BLAST_PW_ALIGN_OUT} "
        "-word_size 3 -matrix blosum62 -gapopen 11 -gapextend 1 "
        "-outfmt 0 -comp_based_stats 0 -seg no -soft_masking false"
