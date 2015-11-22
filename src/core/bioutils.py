"""Genomic sequence processing utilities."""

from base64 import b64encode
from subprocess import Popen
from subprocess import PIPE

from Bio import motifs
from Bio import SeqUtils
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna


def GC(seq):
    """GC content of the sequence"""
    return SeqUtils.GC(seq)


def reverse_complement(seq):
    """Returns the reverse complement of the given sequence."""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement_dict[base] for base in reversed(seq)])


def to_fasta(seqs):
    """Puts given sequences into the FASTA format."""
    fasta_str = ""
    for i, seq in enumerate(seqs):
        fasta_str += ">instance%d\n" % i + seq + "\n"
    return fasta_str


def weblogo(sequences):
    """Generates the sequence logo for the given sequences.

    Uses weblogo program that is locally installed.
    """
    al = to_fasta(sequences)
    weblogo_path = '/usr/local/bin/weblogo'
    p = Popen([weblogo_path, '-F', 'png', '-s', 'LARGE', '-c',
               'classic', '--errorbars', 'YES'],
              stdout=PIPE, stdin=PIPE, stderr=PIPE, close_fds=True)
    stdout_data, stderr_data = p.communicate(input=al)
    return stdout_data


def weblogo_uri(sequences):
    """Generates the weblogo and returns its URI."""
    image_data = weblogo(sequences)
    encoded = b64encode(image_data)
    mime = 'image/png'
    return 'data:' + mime + ';' + 'base64,' + encoded


def build_motif(sites, alphabet=unambiguous_dna):
    """Creates Biopython Motif object from the given sites."""
    instances = [Seq(site, alphabet) for site in sites]
    motif = motifs.create(instances)
    motif.pseudocounts = 0.5
    return motif


def pssm_search(pssm, sequence, threshold=0.0):
    """Finds hits with PSSM score above the threshold."""
    # TODO(sefa): Biopython pssm search doesn't work for ambiguous
    # sequences. Fix this function to handle genome sequences having ambiguous
    # bases.

    dna_sequence = Seq(''.join(b if b in 'ACGT' else 'A' for b in sequence),
                       alphabet=unambiguous_dna)
    scores = pssm.calculate(dna_sequence)
    rc_scores = pssm.reverse_complement().calculate(dna_sequence)
    hits = []
    for pos, (score, rc_score) in enumerate(zip(scores, rc_scores)):
        if score > threshold or rc_score > threshold:
            max_score = max(score, rc_score)
            strand = 1 if score > rc_score else - 1
            hits.append((pos, strand, max_score))
    return hits
                                 
    
