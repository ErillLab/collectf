"""
Given RefSeq accession numbers of the TF-instances in CollecTF, find their
UniProt identifiers.
"""

from collections import Counter
import csv
import json
import pickle
import random
import re

from Bio import Entrez
from tqdm import tqdm
import uniprot

#from base import models

Entrez.email = 'sefa1@umbc.edu'

def fetch_ncbi_protein_record(accession):
    """Fetches the protein record from NCBI, for the given accession number."""
    handle = Entrez.efetch(db='protein', id=accession, retmode='xml')
    records = Entrez.read(handle)
    return records[0]

def parse_wp_accession(rec):
    """Extracts the non-redundant WP accession from a protein record."""
    return rec['GBSeq_contig']

def parse_ncbi_taxonomy_id(rec):
    """Extracts the NCBI taxonomy id from the given protein record."""
    # Not safe, but OK for one-time script.
    source_features, = [features for features in rec['GBSeq_feature-table']
                        if features['GBFeature_key'] == 'source']
    tax_feature, = [feature for feature in source_features['GBFeature_quals']
                    if (feature['GBQualifier_name'] == 'db_xref' and
                        feature['GBQualifier_value'].startswith('taxon:'))]
    tax_id = tax_feature['GBQualifier_value'].replace('taxon:', '')
    return tax_id

def refseq_to_wp(rec):
    """Returns a WP accession for a given NP/YP accession."""
    contig = parse_wp_accession(rec)
    wp_acc = re.match(r'join\((WP_\d+.\d):\d..\d+\)', contig).group(1)
    return wp_acc

def fetch_uniprot_records(uniprot_accessions):
    records = uniprot.retrieve(uniprot_accessions)
    return records.rstrip().split('//\n')
    
def parse_uniprot_tax_id(record):
    """Gets the NCBI taxonomy ID of the given protein."""
    return re.search(r'NCBI_TaxID=(\d+)', record).group(1)

def is_uniprot_record_reviewed(record):
    """Returns true if the UniProt protein record is reviewed."""
    status = re.search(r'ID.*(Reviewed|Unreviewed);', record).group(1)
    assert status in ['Reviewed', 'Unreviewed']
    return status == 'Reviewed'

def is_refseq_accession(accession):
    """Returns true if the given accession number is from RefSeq."""
    prefixes = ['NP', 'YP', 'WP']
    return any(accession.startswith(prefix) for prefix in prefixes)

def is_uniprot_accession(accession):
    """Returns true if the given accession is NOT RefSeq accession number."""
    return not is_refseq_accession(accession)

def filter_same_taxon_accessions(mappings):
    """Filters mappings that have same NCBI taxonomy ID for both RefSeq and
    UniProt accession numbers.
    """
    for refseq_acc in tqdm(mappings):
        if mappings[refseq_acc] and len(mappings[refseq_acc]) > 1:
            print refseq_acc, mappings[refseq_acc]
            refseq_rec = fetch_ncbi_protein_record(refseq_acc)
            refseq_tax_id = extract_ncbi_taxonomy_id(refseq_rec)
            uniprot_accs = [uniprot_acc for uniprot_acc in mappings[refseq_acc]
                            if get_uniprot_tax_id(uniprot_acc) == refseq_tax_id]
            if uniprot_acc:
                mappings[refseq_acc] = uniprot_accs
    return mappings

def get_all_proteins():
    """Returns the list of protein accessions in CollecTF."""
    proteins = [protein.protein_accession
                for protein in models.TFInstance.objects.all()]
    json.dump(proteins, open('scripts/data/all_proteins.json', 'w'))
    return proteins

def refseq_to_uniprot(refseq_accessions):
    """Maps the given RefSeq accession to UniProt ID."""
    return uniprot.map(refseq_accessions, f='P_REFSEQ_AC', t='ACC')

def refseq_to_uniprot_batch(refseq_accs):
    """Given a list of RefSeq accessions, returns the mapping to UniProt. """
    all_mappings = {}
    mapping = refseq_to_uniprot(refseq_accs)
    for refseq_acc in tqdm(refseq_accs):
        all_mappings[refseq_acc] = mapping.get(refseq_acc, None)
        if (not all_mappings[refseq_acc] and
            (refseq_acc.startswith('NP') or refseq_acc.startswith('YP'))):
            wp_acc = refseq_to_wp(refseq_acc)
            all_mappings[refseq_acc] = {wp_acc}
            mapping = refseq_to_uniprot(wp_acc)
            all_mappings[wp_acc] = mapping.get(wp_acc, None)
    return all_mappings

def mock_get_all_proteins():
    proteins =[u'NP_799324', u'NP_231738', u'YP_006516164',
               u'YP_006516595', u'NP_417816']
    return proteins

def migrate_to_uniprot():
    proteins = get_all_proteins()
    all_mappings = refseq_to_uniprot_batch(proteins)
    pickle.dump(all_mappings, open('all_mappings.pkl', 'w'))

def run():
    get_all_proteins()
