"""Genomic sequence processing utilities."""

def reverse_complement(seq):
    """Returns the reverse complement of the given sequence."""
    complement_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([complement_dict[base] for base in reversed(seq)])
