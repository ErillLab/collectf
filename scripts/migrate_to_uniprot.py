"""
Given RefSeq accession numbers of the TF-instances in CollecTF, find their
UniProt identifiers.
"""

from collections import Counter
import csv
import random
import re

from Bio import Entrez
from tqdm import tqdm
import uniprot

#from base import models

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

def extract_ncbi_taxonomy_id(rec):
    """Extracts the NCBI taxonomy id from the given protein record."""
    # Not safe, but OK for one-time script.
    source_features, = [features for features in rec['GBSeq_feature-table']
                        if features['GBFeature_key'] == 'source']
    tax_feature, = [feature for feature in source_features['GBFeature_quals']
                    if (feature['GBQualifier_name'] == 'db_xref' and
                        feature['GBQualifier_value'].startswith('taxon:'))]
    tax_id = tax_feature['GBQualifier_value'].replace('taxon:', '')
    return tax_id

def refseq_to_wp(refseq_accession):
    """Returns a WP accession for a given NP/YP accession."""
    contig = extract_wp_accession(fetch_ncbi_protein_record(refseq_accession))
    wp_acc = re.match(r'join\((WP_\d+.\d):\d..\d+\)', contig).group(1)
    return wp_acc

def get_uniprot_tax_id(uniprot_accession):
    """Gets the NCBI taxonomy ID of the given protein."""
    record = uniprot.retrieve(uniprot_accession)
    return re.search(r'NCBI_TaxID=(\d+)', record).group(1)

class Protein:
    """Class for RefSeq->UniProt accession conversion."""
    def __init__(self, refseq_accession):
        self.refseq_accessions = [refseq_accession]
        self.uniprot_accessions = []
        self.uniprot_accessions_with_same_taxid = []

    def find_uniprot_accessions(self):
        """Finds all mapping UniProt accessions of the protein."""
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
            self.find_uniprot_accessions()

    def find_uniprot_accession_with_same_taxon(self):
        """If there are more than one UniProt mappings, find the one that has
        the same NCBI taxonomy ID with the protein of the RefSeq accession.
        """
        if len(self.uniprot_accessions) > 1:
            refseq_rec = fetch_ncbi_protein_record(self.refseq_accessions[0])
            refseq_tax_id = extract_ncbi_taxonomy_id(refseq_rec)
            #print self.refseq_accessions[0], refseq_tax_id
            for uniprot_acc in self.uniprot_accessions:
                uniprot_tax_id = get_uniprot_tax_id(uniprot_acc)
                #print uniprot_acc, uniprot_tax_id
                if refseq_tax_id == uniprot_tax_id:
                    self.uniprot_accessions_with_same_taxid.append(uniprot_acc)

def get_all_proteins():
    """Returns the list of protein accessions in CollecTF."""
    return [Protein(protein.protein_accession)
            for protein in models.TFInstance.objects.all()]

def run():
    uniprot_mapping_writer = csv.writer(open('refseq2uniprot.csv', 'w'))
    uniprot_mapping_writer.writerow(['RefSeq',
                                     'UniProt (all)',
                                     'UniProt (w same taxid)'])
    acc_sep = '/'               # accession number separator
    proteins = get_all_proteins()
    for protein in tqdm(proteins):
        protein.find_uniprot_accessions()
        protein.find_uniprot_accession_with_same_taxon()
        uniprot_mapping_writer.writerow([
            acc_sep.join(protein.refseq_accessions),
            acc_sep.join(protein.uniprot_accessions),
            acc_sep.join(protein.uniprot_accessions_with_same_taxid)])
                                     

