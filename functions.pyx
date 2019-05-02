# iterate over k-mers in a sequence
def iter_kmers(str seq, int k):
    for i in range(len(seq) - k + 1):
        yield i, seq[i:(i + k)]

# build fuzzy hash map
def build_fuzzy_map(str seq, int k, list seed_exact_idxs,
                    list seed_fuzzy_idxs):
    fuzzy_map = {}
    for kmer_pos, kmer in iter_kmers(seq, k):
        kmer_exact = ''.join([kmer[x] for x in seed_exact_idxs])
        kmer_fuzzy = ''.join([kmer[x] for x in seed_fuzzy_idxs])
        if kmer_exact not in fuzzy_map:
            fuzzy_map[kmer_exact] = {kmer_fuzzy: [kmer_pos]}
        elif kmer_fuzzy not in fuzzy_map[kmer_exact]:
            fuzzy_map[kmer_exact][kmer_fuzzy] = [kmer_pos]
        else:
            fuzzy_map[kmer_exact][kmer_fuzzy].append(kmer_pos)
    return fuzzy_map
