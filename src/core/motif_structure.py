from Bio import motifs
import itertools

def build_motif(sites):
    """Builds a Biopython motifs.Motif object out of given sites."""
    motif = motifs.create(sites)
    motif.pseudocounts = 0.8   # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2647310/
    return motif


def slice_sites(sites, start, end):
    """Slices each site."""
    return [site[start:end] for site in sites]


def score_site(pssm, site):
    """Scores the given site with the given PSSM."""
    return sum(pssm[site[i]][i] for i in range(len(site)))


def score_sites(pssm, sites):
    """Computes the average PSSM score of a list of sites."""
    return sum(score_site(pssm, site) for site in sites) / len(sites)


def direct_repeat(seq):
    """Match for a sequence in a direct-repeat pattern."""
    return seq


def inverted_repeat(seq):
    """Match for a sequence in an inverted-repeat pattern."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement[b] for b in seq[::-1])


def find_pattern(sites, self_score_ratio_threshold=0.6,
                 kmer_pair_score_ratio_threshold=0.3):
    """Finds pattern in a motif."""
    k = 4
    sites_len = len(sites[0])
    all_kmers = [{'start': i, 'end': i+k, 'seqs': slice_sites(sites, i, i+k)}
                 for i in range(sites_len-k+1)]
    # Compute self-PSSM scores
    for kmer in all_kmers:
        motif = build_motif(kmer['seqs'])
        kmer['pssm'] = motif.pssm
        kmer['self_score'] = score_sites(kmer['pssm'], kmer['seqs'])

    max_self_score = max(kmer['self_score'] for kmer in all_kmers)
    all_kmers = [kmer for kmer in all_kmers
                 if kmer['self_score'] > self_score_ratio_threshold*max_self_score]


    pattern = (0, 'Single box')
    for kmer_a, kmer_b in itertools.combinations(all_kmers, 2):
        if not (kmer_a['start'] >= kmer_b['end'] or
                kmer_b['start'] >= kmer_a['end']):
            continue

        if kmer_a['self_score'] < kmer_b['self_score']:
            kmer_a, kmer_b = kmer_b, kmer_a

        # Look for direct-repeat
        score = score_sites(kmer_a['pssm'], kmer_b['seqs'])
        if (score > kmer_pair_score_ratio_threshold*kmer_a['self_score'] and
            score > pattern[0]):
            pattern = (score, 'Direct repeat', kmer_a['start'], kmer_b['start'])

        # Look for inverted-repeat
        score = score_sites(
            kmer_a['pssm'], [inverted_repeat(site) for site in kmer_b['seqs']])
        if (score > kmer_pair_score_ratio_threshold*kmer_a['self_score'] and
            score > pattern[0]):
            pattern = (score, 'Inverted repeat', kmer_a['start'], kmer_b['start'])

    return pattern[1]
