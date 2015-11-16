"""Functions to fetch information from NCBI via Entrez."""

import urllib

from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'sefa1@umbc.edu'


class EntrezException(Exception):
    pass


def get_pubmed(pmid):
    """Retrieves the PubMed publication from NCBI database."""
    try:
        handle = Entrez.esummary(db='pubmed', id=pmid)
        record = Entrez.read(handle)
        return record[0]
    except RuntimeError:
        raise EntrezException


def get_genome(accession):
    """Retrieves the genome from RefSeq database."""
    try:
        handle = Entrez.efetch(db='nuccore', id=accession,
                               retmode='gbwithparts', rettype='text')
        record = SeqIO.read(handle, 'gb')
        handle.close()
        return record
    except urllib.error.HTTPError:
        raise EntrezException
