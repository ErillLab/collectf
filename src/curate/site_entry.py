"""This file contains functions that are used for parsing and processing the
text entered on site-entry step.

First, the text is parsed to find out whether it is sequence-based or
coordinate-based. In addition, the text will be checked to see if it contains
any quantitative information. All fields must be separated by space or tab.  """

from django.utils.safestring import mark_safe

import re
import regex
import StringIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm

from core import bioutils
from core import misc


class Match:
    """Match class for site matches on the genome"""
    def __init__(self, genome, seq, reported_seq, start, end, strand):
        """Initializer for the match object.
        - genome: The genome object that match is in.
        - seq: The sequence in the genome
        - reported_seq: The original sequence reported in the paper. It can be
          different from the matched seq
        - start,end,strand: location of sequence in the genome
        """
        self.genome = genome  # the genome
        self.seq = str(seq)  # the sequence in the genome
        self.reported_seq = str(reported_seq)  # the sequence as reported
        self.start = start
        self.end = end
        self.strand = strand
        self.regulated_genes = []

    def set_nearby_genes(self, dist_th=150):
        """Given a match, find the genes nearby.
        The method to find nearby genes is as follows:
        Find the genes that are left (and right) of the site. Starting from the
        closest gene, scan the genome to the left (and right), until no gene is
        find in the next <dist_th> base-pair."""
        genes = self.genome.get_genes()
        # make sure all genes are sorted by start position
        assert all(ga.start <= gb.start for ga, gb in zip(genes, genes[1:]))
        left_genes = filter(lambda g: g.end < self.start, genes)
        right_genes = filter(lambda g: g.start > self.end, genes)
        # genes overlap with site
        nearby_genes = genes[len(left_genes):len(genes)-len(right_genes)]
        # left
        if left_genes:
            ngi, nearby_gene = min(enumerate(left_genes),
                                   key=lambda x: dist(x[1], self))
            if nearby_gene.strand == -1:
                nearby_genes.append(nearby_gene)
                while (ngi > 0 and
                       left_genes[ngi-1].strand == -1 and
                       dist(left_genes[ngi-1], left_genes[ngi]) < dist_th):
                    nearby_genes.append(left_genes[ngi-1])
                    ngi -= 1
        # right
        if right_genes:
            ngi, nearby_gene = min(enumerate(right_genes),
                                   key=lambda x: dist(x[1], self))
            if nearby_gene.strand == 1:
                nearby_genes.append(nearby_gene)
                while (ngi < len(right_genes)-1 and
                       right_genes[ngi+1].strand == 1 and
                       dist(right_genes[ngi], right_genes[ngi+1]) < dist_th):
                    nearby_genes.append(right_genes[ngi+1])
                    ngi += 1

        if not nearby_genes:
            # if there is no site nearby, just add the nearest one
            nearby_genes.append(min(genes, key=lambda x: dist(x, self)))
        self.nearby_genes = nearby_genes

    def clear_regulated_genes(self):
        """Clear regulated genes by TF binding to the site"""
        self.regulated_genes = []

    def set_regulated_genes(self, gs):
        """Add regulated genes."""
        self.regulated_genes = gs

    def is_exact(self):
        return self.seq == self.reported_seq

    def pprint(self, diagram_hover=True):
        """Given a match object, make the HTML snippet to display it
        properly."""
        # Make the diagram and gene function table
        extra = (self.match_diagram() +
                 """<table class="table table-condensed small">
                 <thead><tr><th>locus tag</th><th>gene
                 name</th><th>function</th></tr></thead> <tbody>""")
        for g in self.nearby_genes:
            extra += ("<tr><td>%s</td><td>%s</td><td>%s</td></tr>" %
                      (g.locus_tag, g.name, g.description))
        extra += ('</tbody>' + "</table><br/>")
        return_str = ""
        if self.is_exact():
            return_str += ('<div class="sequence">%s<br/> %s[%d,%d] %s</div>' %
                           (self.seq, '+' if self.strand == 1 else '-',
                            self.start+1, self.end+1,
                            self.genome.genome_accession))
        else:
            return_str += self.print_alignment(self.reported_seq, self.seq)

        if diagram_hover:
            return mark_safe(return_str + extra)

        return mark_safe(return_str)

    def match_diagram(self):
        """This method is called during curation submission. When the reported
        sites are given, they are searched in the genome and matched are
        displayed. For display, site match diagram is created and presented."""
        gdd = GenomeDiagram.Diagram("Site match diagram")
        gdt_features = gdd.new_track(1, greytrack=False)
        gds_features = gdt_features.new_set()
        # draw genes
        for gene in self.nearby_genes:
            feature = SeqFeature(FeatureLocation(gene.start+1, gene.end+1),
                                 strand=gene.strand)
            gds_features.add_feature(
                feature,
                name=gene.name,
                label=True,
                label_size=10,
                label_angle=0 if gene.strand == 1 else 180,
                label_position='middle',
                sigil='ARROW',
                arrowshaft_height=1.0,
                color=colors.lightblue)
        # Draw binding site
        feature = SeqFeature(FeatureLocation(self.start+1, self.end+1),
                             strand=self.strand)
        gds_features.add_feature(feature, color=colors.red, name='site',
                                 label=False, label_size=12)
        gdd.draw(format='linear', fragments=1,
                 start=min(map(lambda g: g.start+1, self.nearby_genes))-150,
                 end=max(map(lambda g: g.end+1, self.nearby_genes))+150,
                 pagesize=(2*cm, 12*cm))
        return mark_safe('<div>' + gdd.write_to_string('svg') + '</div>')

    def print_alignment(self, seqa, seqb):
        """Given two sequences, pairwise align them and output HTML for curation
        exact/inexact site match steps"""
        assert len(seqa) == len(seqb)
        s = []
        s.append('<span class="sequence">')
        s.append("%s<br/>" % seqa)
        s.append(''.join('|' if seqa[i]==seqb[i] else '&nbsp;' for i in range(len(seqa))) + '<br/>')
        s.append('%s %s[%d,%d] (%s)</br>' % (seqb, '+' if self.strand==1 else '-',
                                             self.start, self.end, self.genome.genome_accession))
        s.append('</span>')
        return mark_safe( ''.join(s))

    def __repr__(self):
        return "%s %s (exact? %s)" % (self.seq, self.reported_seq, self.is_exact())


class Site:
    @property
    def key(self):
        return self._key

    def set_nearby_genes_for_all_matches(self):
        """Find nearby genes for all matches for a site."""
        for match in self.exact_matches:
            match.set_nearby_genes()
        for match in self.soft_matches:
            match.set_nearby_genes()

    def populate_match_choices(self, add_no_valid_opt, match_type):
        """For a given site and its all matches, populate Django field choices.
        If match_type is 'exact_only', then populate choices for exact matches
        only. Similarly, if the match_type is 'inexact_only', then populate
        choices for soft-search matches only."""
        if match_type == 'exact_only':
            matches = self.exact_matches
        elif match_type == 'inexact_only':
            matches = self.soft_matches
        else:
            assert False, "invalid match_type option"

        choices = [(i, match.pprint()) for (i, match) in enumerate(matches)]
        if add_no_valid_opt:
            choices.append((None, "No valid match."))
        return choices

    def get_exact_matches(self):
        """Return all exact matches"""
        return self.exact_matches

    def set_exact_match(self, match_id):
        """Given a match id (one of site's possible exact matches), match the
        genome location to the binding site reported in the paper"""
        self.matched = self.exact_matches[int(match_id)]

    def set_soft_match(self, match_id):
        """Given a match id (one of site's possible soft matches), match the
        genome location to the binding site reported in the paper"""
        self.matched = self.soft_matches[int(match_id)]

    def get_match(self):
        """Return the matched location (Match object)"""
        assert self.matched, "Not matched yet."
        return self.matched

    def is_matched(self):
        return bool(self.matched)

    def set_qval(self, qval):
        """Set quantitative value."""
        self.qval = qval

    def set_TF_function(self, TF_function):
        """Set TF function for the site."""
        self.TF_function = TF_function

    def set_TF_type(self, TF_type):
        """Set TF type for the site."""
        self.TF_type = TF_type

    def clear_techniques(self):
        """Clear experimental techniques used to determine the site."""
        self.techniques = []

    def add_technique(self, t):
        """Add an experimental technique to the set of techniques used to
        determine the site."""
        self.techniques.append(t)

    def match_peak_data(self, peaks):
        """If the curation mode is high-throughput, the list of peaks are read
        from the 'peaks' field. In case of the quantitative values that are
        reported in this field, each site in the 'sites' field is attempted to
        match one of the peaks in the 'peaks' field to link the associated
        quantitative value"""
        print 'match_peak_data'
        # If peaks can not have quantitative value, just skip this function
        if not peaks[0].qval:
            return
        # If not matched to any site in the genome, the site will be saved as
        # not-annotated. In this case, the quantitative value is not necessary,
        # so skip it.
        if not self.matched:
            return

        match = self.get_match()
        for peak in peaks:
            if hasattr(peak, 'start'): # coordinate-based peaks
                if (min(peak.start, peak.end) <= match.start and
                    max(peak.start, peak.end) >= match.end):
                    self.qval = peak.qval
                    break
            else: # sequence-based peaks
                if (match.seq in peak.seq or
                    bioutils.reverse_complement(match.seq) in peak.seq):
                    self.qval = peak.qval
                    break

class SequenceSite(Site):
    """Class definition for sites that are initialized with the sequence"""
    def __init__(self, id, seq, qval=None):
        self._key = id
        if any(nuc not in 'ACGT' for nuc in seq):
            raise
        self.seq = seq
        self.qval = qval
        self.exact_matches = []
        self.soft_matches = []
        self.matched = None

    def __repr__(self):
        return "%s [%.2f]" % (self.seq, self.qval if self.qval else 0)

    def pprint(self):
        return mark_safe('<span class="sequence">%s</span>' % self.seq)

    def search_exact_match(self, genomes):
        """Search the genome and find exact matches on the genome."""
        self.exact_matches = []
        self.soft_matches = []
        self.matched = None
        for genome in genomes:
            self.exact_matches.extend(self.locate_seq(genome))
        # find nearby genes for all matches
        self.set_nearby_genes_for_all_matches()

    def locate_seq(self, genome):
        """Search the site sequence in the genome, both strands"""
        matches = (self._locate_seq_strand(genome, 1) +
                   self._locate_seq_strand(genome, -1))
        return matches

    def _locate_seq_strand(self, genome, strand):
        """Search the site sequence only in one strand of the genome."""
        search_seq = self.seq if strand==1 else bioutils.reverse_complement(self.seq)
        matches = []
        i = genome.sequence.find(search_seq)
        while i >= 0:
            matches.append(Match(
                genome, self.seq, self.seq, i, i+len(search_seq)-1, strand))
            i = genome.sequence.find(search_seq, i+1)
        return matches

    def search_soft_match(self, genomes, motif=None):
        """Search the genome and find in-exact matches on the genome."""
        # overwrite exact_matches (if any) with soft-search results
        self.soft_matches = []
        self.matched = None
        for genome in genomes:
            self.soft_matches.extend(self.soft_locate_seq(genome))
        # Find nearby genes for all soft-matches.
        self.set_nearby_genes_for_all_matches()

    def soft_locate_seq(self, genome, mismatch_th=2, motif=None):
        """(1) Find all sequences in the genome such that they are at most
        <mismatch_th> far.  (2) Run PSSM built with sites in the DB and sort
        search results by score."""
        matches = (self._soft_locate_seq_strand(genome, 1, mismatch_th) +
                   self._soft_locate_seq_strand(genome, -1, mismatch_th))
        return matches

    def _soft_locate_seq_strand(self, genome, strand, mismatch_th):
        """Soft search on one strand only. Return all matches that are upto
        <mismatch_th> away from the reported sequence."""
        # First, find all sequences on the genome that are similar enough
        genome_sequence = genome.sequence
        search_seq = self.seq if strand == 1 else bioutils.reverse_complement(self.seq)
        pattern = regex.compile('(%s){1<=s<=%d}' % (search_seq, mismatch_th))
        matches = []
        for m in regex.finditer(pattern, genome_sequence, overlapped=True):
            matched_seq = m.group() if strand==1 else bioutils.reverse_complement(m.group())
            match = Match(genome, matched_seq, self.seq, m.span()[0], m.span()[1]-1, strand)
            matches.append(match)
        return matches

    def __repr__(self):
        return "%s" % self.seq

class CoordinateSite(Site):
    """Class definition for sites that are initialized using coordinates."""
    def __init__(self, id, start, end, qval=None):
        self._key = id
        self.start = start
        self.end = end
        self.qval = qval
        self.exact_matches = []
        self.soft_matches = []
        self.matched = None

    def search_exact_match(self, genomes):
        """It performs the same job with SequenceSite search_exact_match
        function, which is finding the exact matches in the genome. Since
        coordinates are given here, all this function does is to get the
        corresponding region from the genome."""
        self.exact_matches = []
        self.soft_matches = []
        self.matched = None
        start = self.start
        end = self.end
        for genome in genomes:
            genome_sequence = genome.sequence
            if self.start <= self.end:
                # there is only one match
                self.seq = genome_sequence[start:end+1]
                self.exact_matches = [Match(genome, self.seq, self.seq, start, end, 1)]
            else:
                self.seq = bioutils.reverse_complement(genome_sequence[end:start+1])
                self.exact_matches = [Match(genome, self.seq, self.seq, end, start, -1)]

        # For all matches, get nearby genes
        self.set_nearby_genes_for_all_matches()

def parse_fasta(text):
    """Parse fasta file and return list of sites"""
    l = SeqIO.parse(StringIO.StringIO(text), 'fasta')
    seqs = [item.seq.tostring() for item in l]
    return [SequenceSite(i,seq) for i,seq in enumerate(seqs)]

def parse_seq(text):
    """Parse text that contains a list of sequences, one per line"""
    seqs = [l.strip() for l in text.split('\n') if l]
    return [SequenceSite(i,seq) for i,seq in enumerate(seqs)]

def parse_seq_with_qval(text):
    """Parse text that contains a list of sequences and associated quantitative
    values. There should be one pair of site and quantitative value on each line
    of text."""
    lines = re.split('[\r\n]+', text)
    seqs = [line.split()[0] for line in lines]
    quantitative_values = map(float, [line.split()[1] for line in lines])
    assert len(seqs) == len(quantitative_values)
    return [SequenceSite(i, seq, qval) for i,(seq,qval) in enumerate(zip(seqs, quantitative_values))]

def parse_coords(text):
    """Parse text that contains list of coordinates, one per line. Each line
    must contain a pair of numbers, denoting start and end positions for the
    site, respectively."""
    coordinates = [re.split('[\t ]+', line)
                   for line in re.split('[\r\n]+', text)]
    return [CoordinateSite(i, int(coord[0])-1, int(coord[1])-1)
            for i, coord in enumerate(coordinates)]

def parse_coords_with_qval(text):
    """Parse text that contains list of coordinates, one per line. Additionally,
    each line has a quantitative value associated with the coordinates"""
    coordinates = [re.split('[\t ]+', line)
                   for line in re.split('[\r\n]+', text)]
    return [CoordinateSite(i, int(coord[0])-1, int(coord[1])-1, float(coord[2]))
            for i, coord in enumerate(coordinates)]

def parse_input(text):
    """Parse text of reported sites.
    It can contain sequences or coordinates, optionally with quantitative
    values.
    """
    sites = None
    if text.strip().startswith('>'): # fasta format
        sites = parse_fasta(text)
    else: # it can be in sequences/coordinates (and quantitative values)
        site_lines = re.split('[\r\n]+', text)
        if len(site_lines[0].split()) == 3: # should be coordinates with q values
            sites = parse_coords_with_qval(text)
        elif len(site_lines[0].split()) == 2:
            # Can be either sequence with q values or
            # coordinates without q values
            a,b = site_lines[0].split()
            if a.isalpha() and misc.is_float(b):
                sites = parse_seq_with_qval(text)
            elif a.isdigit() and b.isdigit():
                sites = parse_coords(text)
        elif len(site_lines[0].split()) == 1:
            # it should be only sequences
            sites = parse_seq(text)
    return sites

def dist(a, b):
    """Given two genes (or gene and match), return distance"""
    return max(a.start, b.start) - min(a.end, b.end)
