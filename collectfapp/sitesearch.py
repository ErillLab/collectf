from Bio import SeqIO
from Bio.Seq import Seq
import StringIO
from collections import namedtuple
import regex

# Some namedtuple declarations
Match = namedtuple('Match', 'seq start end strand')
SiteMatch = namedtuple('SiteMatch', 'match nearby_genes')
def print_match(m):
    strand = '+' if m.strand == 1 else '-'
    return u'%s, %s[%d, %d]' % (m.seq, strand, m.start, m.end)

def parse_site_input(text):
    """Parse text of reported sites. It can be either in FASTA format, or plain
    list of site sequences"""
    if text.strip().startswith('#'): # FASTA format
        l = bio.SeqIO.parse(StringIO(text), 'fasta')
        sites = [item.seq.tostring() for item in l]
    else: # plain '\n' separated list of sequences
        sites = [l for l in text.split() if l]
    return dict(enumerate(sites)) # list of (sid, site)

def locate_site_strand(genome_seq, site_seq, strand):
    """Site search on genome sequence"""
    matches = []
    i = genome_seq.find(site_seq)
    while i >= 0:
        m = Match(seq=site_seq, start=i, end=i+len(site_seq)-1, strand=strand)
        matches.append(m)
        i = genome_seq.find(site_seq, i+1)
    return matches

def locate_site(genome_seq, site_seq):
    # find exact matches of site on the genome sequence
    # search both strands
    reverse_site_seq = Seq(site_seq).reverse_complement().tostring()
    matches = locate_site_strand(genome_seq, site_seq, strand=1) + \
              locate_site_strand(genome_seq, reverse_site_seq, strand=-1)
    return matches

def soft_locate_site_strand(genome_seq, site_seq, strand, pattern):
    """Soft seatch on one strand only"""
    # first find all sequences on genome that are similar enough
    matches = []
    for m in regex.finditer(pattern, genome_seq, overlapped=True):
        match = Match(seq=m.group(), start=m.span()[0], end=m.span()[1], strand=1)
        matches.append(match)
    return matches

def soft_locate_site(genome_seq, site_seq, mismatch_th=2):
    """SOFT search for site in genome sequence.
    Since, PSSM search over whole genome is _very_ expensive,
    - find all sequences on the genome such that they are at most mismatch_th far.
    - run PSSM over those sequences and sort them by score.
    In future Pat's suffix array module may be used."""
    reverse_site_seq = Seq(site_seq).reverse_complement().tostring()
    # allow upto <mismatch_th> substitutions
    p = regex.compile('(%s){1<=s<=%d}' % (site_seq, mismatch_th))
    matches = soft_locate_site_strand(genome_seq, site_seq, strand=1, pattern=p) + \
              soft_locate_site_strand(genome_seq, reverse_site_seq, strand=-1, pattern=p)
    return matches
    
def dist(a, b):
    """Given two genes (or gene and Match), return distance"""
    return max(a.start, b.start) - min(a.end, b.end)

def locate_nearby_genes(genes, site_loc, dist_th=50):
    """Given the list of models.Gene objects, locate genes close to the
    site_loc. Return list of nearby genes."""
    # make sure all genes are sorted by start position
    assert all(ga.start <= gb.start for ga,gb in zip(genes, genes[1:]))

    nearby_genes = []
    # Instead of checking every gene, perform binary search to find nearby ones.
    low, high = 0, len(genes)
    intersect = False  # one of genes and site_loc intersects?
    while low < high and not intersect:
        mid = (low + high) // 2
        if site_loc.start > genes[mid].end:
            low = mid + 1
        elif site_loc.end < genes[mid].start:
            high = mid
        else: # it means genes[mid] intersects with site_loc
            intersect = True
            nearby_genes.append(genes[mid])
    # For operons, move left/right while two genes are close enough
    lhs_index, rhs_index = (mid-1, mid+1) if intersect else (low-1, low)
    # to the left
    if lhs_index >= 0:
        nearby_genes.append(genes[lhs_index])
    while lhs_index>0 and dist(genes[lhs_index-1], genes[lhs_index]) < dist_th:
        nearby_genes.append(genes[lhs_index-1])
        lhs_index -= 1
    # to the right
    if rhs_index < len(genes):
        nearby_genes.append(genes[rhs_index])
    while (rhs_index < len(genes)-1 and
           dist(genes[rhs_index], genes[rhs_index+1]) < dist_th):
        nearby_genes.append(genes[rhs_index+1])
        rhs_index += 1
    return nearby_genes    

def match_all(genome, genes, sites, exact):
    """For all sites in list, find exact/soft matches on the genome"""
    # For each site, store a SiteMatch namedtuple object.
    match_function = locate_site if exact else soft_locate_site
    site_matches = {} # dictionary of {sid: {mid: SiteMatch}}
    for sid, site in sites.items():
        site_matches[sid] = {} # initialize
        matches = match_function(genome.sequence, site)
        # for each match, find nearby genes
        for mid, match in enumerate(matches):
            nearby_genes = locate_nearby_genes(genes, match)
            site_matches[sid][mid] = SiteMatch(match=match, nearby_genes=nearby_genes)
    return site_matches
            
        
