"""Genomic sequence processing utilities."""

from base64 import b64encode
from subprocess import Popen
from subprocess import PIPE

from Bio import motifs
from Bio.Seq import Seq

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


def build_motif(sites):
    """Creates Biopython Motif object from the given sites."""
    instances = [Seq(site) for site in sites]
    return motifs.create(instances)

