# Given RefSeq accession numbers of the TF-instances in CollecTF, find their
# UniProt identifiers.

import re
from Queue import Queue

import uniprot
from Bio import Entrez
from tqdm import tqdm

from base import models

Entrez.email = 'sefa1@umbc.edu'

def refseq_to_uniprot(refseq_accession):
    """Maps the given RefSeq accession to UniProt ID."""
    return uniprot.map(refseq_accession, f='P_REFSEQ_AC', t='ACC')

def fetch_ncbi_protein_record(accession):
    """Fetches the protein record from NCBI, for the given accession number."""
    handle = Entrez.efetch(db='protein', id=accession, retmode='xml')
    records = Entrez.read(handle)
    return records[0]

def extract_wp_accession(rec):
    """Extracts the non-redundant WP accession from a protein record."""
    return rec['GBSeq_contig']

def refseq_to_wp(refseq_accession):
    """Returns a WP accession for a given NP/YP accession."""
    contig = extract_wp_accession(fetch_ncbi_protein_record(refseq_accession))
    wp_acc = re.match(r'join\((WP_\d+.\d):\d..\d+\)', contig).group(1)
    return wp_acc

class Protein:
    def __init__(self, refseq_accession):
        self.refseq_accessions = [refseq_accession]
        self.uniprot_accessions = []

    def find_uniprot_accession(self):
        """Finds the UniProt accession of the protein."""
        refseq_acc = self.refseq_accessions[-1]
        # Try to gets its UniProt mapping
        mapping = refseq_to_uniprot(refseq_acc)
        if refseq_acc in mapping:
            self.uniprot_accessions.extend(mapping[refseq_acc])
        else:
            if refseq_acc.startswith('WP'):
                return
            # Convert NP/YP record to its WP counterpart and try accession
            # number with WP for UniProt mapping.
            wp_refseq_acc = refseq_to_wp(refseq_acc)
            self.refseq_accessions.append(wp_refseq_acc)
            self.find_uniprot_accession()
            
def get_all_proteins():
    """Returns the list of protein accessions in CollecTF."""
    return [Protein(protein.protein_accession)
            for protein in models.TFInstance.objects.all()]

def run():
    proteins = get_all_proteins()
    for protein in proteins:
        protein.find_uniprot_accession()
        print protein.refseq_accessions, protein.uniprot_accessions
