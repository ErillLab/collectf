from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Motif
from Bio.Alphabet import IUPAC
import StringIO
from collections import namedtuple
import regex
from django.utils.safestring import mark_safe
from templatetags import utils

from baseapp.templatetags import gene_diagram


# Some namedtuple declarations
Match = namedtuple('Match', 'seq start end strand')
SiteMatch = namedtuple('SiteMatch', 'match nearby_genes')

def print_alignment(seqa, seqb):
    """Given two sequences, pairwise align them and output HTML for curation
    exact/inexact site match steps"""
    assert len(seqa) == len(seqb)
    s = []
    s.append('<span class="sequence">')
    s.append("%s<br/>" % seqa)
    s.append(''.join('|' if seqa[i]==seqb[i] else '&nbsp;' for i in range(len(seqa))) + '<br/>')
    s.append('%s</br>' % seqb)
    s.append('</span>')
    return mark_safe( ''.join(s))

def print_site_match(reported_site, m, is_exact):
    """Given a match object, make the html snippet to display it nice.  I am not sure
    this is a proper solution (HTML in python), couldn't find a better&easier way to
    do though"""
    
    strand = '+' if m.match.strand==1 else '-'
    nearby_genes = [g.locus_tag + (' (%s)' % g.name if g.name != g.locus_tag else '')
                    for g in m.nearby_genes]
    s = ""
    if is_exact:
        s += ('<span class="sequence"> %s %s(%d, %d)</span><br/>' %
              (m.match.seq, '+' if m.match.strand == 1 else '-',  m.match.start, m.match.end))
    else:
        s += print_alignment(reported_site, m.match.seq)

    s += (gene_diagram.site_match_diagram(m) +
          '<table class="table table-condensed">' +
          '<thead><tr><th>locus tag</th><th>gene name</th><th>function</th></tr></thead>' +
          '<tbody>'
          )
    for g in m.nearby_genes:
        s += "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (g.locus_tag, g.name, g.description)
    s += ('</tbody>' + "</table><br/>")

    return mark_safe(s)  # render newline correctly

def parse_site_input(text):
    """Parse text of reported sites. It can be either in FASTA format, or plain
    list of site sequences"""
    if text.strip().startswith('>'): # FASTA format
        l = SeqIO.parse(StringIO.StringIO(text), 'fasta')
        sites = [item.seq.tostring() for item in l]
    else: # plain '\n' separated list of sequences
        sites = [l for l in text.split() if l]

    # make uppercase
    sites = [site.upper() for site in sites]
    # if any site have ambiguous nucs, return None
    if any(nuc not in "ACGTacgt" for site in sites for nuc in site):
        return None
    return dict(enumerate(sites)) # list of (sid, site)

def reverse_complement(seq):
    return Seq(seq).reverse_complement().tostring()

def locate_site_strand(genome_seq, site_seq, strand):
    """Site search on genome sequence"""
    search_seq = site_seq if strand==1 else reverse_complement(site_seq)
    matches = []
    i = genome_seq.find(search_seq)
    while i >= 0:
        m = Match(seq=site_seq, start=i, end=i+len(site_seq)-1, strand=strand)
        matches.append(m)
        i = genome_seq.find(search_seq, i+1)
    return matches

def locate_site(genome_seq, site_seq):
    # find exact matches of site on the genome sequence
    # search both strands
    matches = locate_site_strand(genome_seq, site_seq, strand=1) + \
              locate_site_strand(genome_seq, site_seq, strand=-1)
    return matches

def soft_locate_site_strand(genome_seq, strand, site_seq, mismatch_th):
    """Soft seatch on one strand only"""
    # first find all sequences on genome that are similar enough
    search_seq = site_seq if strand==1 else reverse_complement(site_seq)
    pattern = regex.compile('(%s){1<=s<=%d}' % (search_seq, mismatch_th))
    matches = []
    # find all matches upto <mismatch_th> substitutions
    for m in regex.finditer(pattern, genome_seq, overlapped=True):
        matched_seq = m.group() if strand==1 else reverse_complement(m.group())
        match = Match(seq=matched_seq, start=m.span()[0], end=m.span()[1]-1, strand=strand)
        matches.append(match)
    return matches

def soft_locate_site(genome_seq, site_seq, mismatch_th=2, motif=None):
    """SOFT search for site in genome sequence.
    Since, PSSM search over whole genome is _very_ expensive,
    - find all sequences on the genome such that they are at most mismatch_th far.
    - run PSSM over those sequences and sort them by score.
    In future Pat's suffix array module may be used."""
    matches = soft_locate_site_strand(genome_seq, 1, site_seq, mismatch_th) + \
              soft_locate_site_strand(genome_seq, -1, site_seq, mismatch_th)
    # sort matches based on PWM
    if not motif:
        return matches

    matches.sort(key=lambda x: score_match(motif, x.seq), reverse=True)
    
    return matches

def build_motif(seqs):
    """Create motif from sequences"""
    m = Motif.Motif(alphabet=IUPAC.unambiguous_dna)
    for seq in seqs:
        m.add_instance(Seq(seq, m.alphabet))
    return m

def score_match(motif, sequence):
    """Given Biopython motif object and a sequence match, return PWM score of the match"""
    return motif.score_hit(sequence, position=0)
    
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

def match_all_exact_coordinates(genome, genes, coordinates):
    """Return sites and matches for chip-seq data. Session Data structures don't make
    much sense, but it is this way to be in the same format with motif-associated
    sites"""
    
    sites = {}
    site_matches = {} # for all coordinates, there will be one exact match
    peak_intensities = {} 
    for cid, coor in enumerate(coordinates):
        start = int(coor[0])
        end = int(coor[1])
        print start,end
        match = Match(seq=genome.sequence[start-1:end-1], start=start, end=end, strand=1)
        nearby_genes = locate_nearby_genes(genes, match)
        sites[cid] = match.seq
        site_matches[cid] = SiteMatch(match=match, nearby_genes=nearby_genes)
        if len(coor) == 3:
            peak_intensities[cid] = float(coor[2])
    return sites, site_matches, peak_intensities
    

def match_all_exact(genome, genes, sites):
    """For all sites in list, find exact matches on the genome"""
    # For each site, store a SiteMatch namedtuple object.
    site_matches = {} # dictionary of {sid: {mid: SiteMatch}}
    for sid, site in sites.items():
        site_matches[sid] = {} # initialize
        matches = locate_site(genome.sequence, site)
        # for each match, find nearby genes
        for mid, match in enumerate(matches):
            nearby_genes = locate_nearby_genes(genes, match)
            site_matches[sid][mid] = SiteMatch(match=match, nearby_genes=nearby_genes)        
    return site_matches

def match_all_soft(genome, genes, sites, exact_sites):
    """For all sites in the list, find soft matches on the genome. Sort set of
    matches of every site."""
    # create motif from exact sites
    motif = build_motif(exact_sites)
    site_matches = {} # dictionary of {sid: {mid: SiteMatch}}
    for sid, site in sites.items():
        site_matches[sid] = {}  # init
        matches = soft_locate_site(genome.sequence, site, mismatch_th=2, motif=motif)
        # for each match, find nearby genes
        for mid, match in enumerate(matches):
            nearby_genes = locate_nearby_genes(genes, match)
            site_matches[sid][mid] = SiteMatch(match=match, nearby_genes=nearby_genes)
    return site_matches
