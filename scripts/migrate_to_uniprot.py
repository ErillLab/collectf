"""
Given RefSeq accession numbers of the TF-instances in CollecTF, find their
UniProt identifiers.
"""

from collections import Counter
import csv
import json
import os
import pickle
import random
import re
from pprint import pprint

from Bio import Entrez
from tqdm import tqdm
from tqdm import trange
import uniprot

#from base import models

#DATA_DIR = '/home/sefa/Dropbox/collectf/scripts/data'
DATA_DIR = '/Users/sefa/Dropbox/collectf/scripts/data'
Entrez.email = 'sefa1@umbc.edu'

def fetch_ncbi_protein_record(accession):
    """Fetches the protein record from NCBI, for the given accession number."""
    handle = Entrez.efetch(db='protein', id=accession, retmode='xml')
    records = Entrez.read(handle)
    return records[0]

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

def parse_wp_accession(rec):
    """Returns a WP accession for a given NP/YP accession."""
    contig = rec['GBSeq_contig']
    wp_acc = re.match(r'join\((WP_\d+.\d):\d..\d+\)', contig).group(1)
    return wp_acc

def fetch_uniprot_records(uniprot_accessions):
    """Returns the list of UniProt records for the given accession numbers."""
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

def is_uniparc_accession(accession):
    """Returns true if the given accession is a UniParc accession number."""
    return accession.startswith('UPI')

def get_all_proteins():
    """Returns the list of protein accessions in CollecTF."""
    json_dump_file = os.path.join(DATA_DIR, 'proteins.json')
    if not os.path.isfile(json_dump_file):
        proteins = [protein.protein_accession
                    for protein in models.TFInstance.objects.all()]
        json.dump(proteins, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def fetch_all_ncbi_protein_records():
    """For all proteins, retrieves the NCBI records."""
    refseq_accessions = get_all_proteins()
    json_dump_file = os.path.join(DATA_DIR, 'ncbi_records.json')
    if not os.path.isfile(json_dump_file):
        all_records = {}
        for refseq in tqdm(refseq_accessions):
            all_records[refseq] = fetch_ncbi_protein_record(refseq)
        json.dump(all_records, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def fetch_all_uniprot_records():
    """For all UniProt accession numbers, retrieves the UniProt records"""
    json_dump_file = os.path.join(DATA_DIR, 'uniprot_records.json')
    if not os.path.isfile(json_dump_file):
        uniprot_mappings = map_refseq_to_uniprot()
        uniprot_accs = [acc for accs in uniprot_mappings.values()
                        for acc in accs if not is_uniparc_accession(acc)]
                        
        all_records = {}
        print "Number of UniProt accession numbers:", len(uniprot_accs)
        
        # Split queries into smaller chunks
        batch_size = 25
        for i in trange(0, len(uniprot_accs), batch_size):
            query_accs = uniprot_accs[i:i+batch_size]
            records = fetch_uniprot_records(query_accs)
            for acc, record in zip(query_accs, records):
                all_records[acc] = record

        json.dump(all_records, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def refseq_to_uniprot(refseq_accessions):
    """Maps the given RefSeq accession to UniProt ID."""
    return uniprot.map(refseq_accessions, f='P_REFSEQ_AC', t='ACC')

def map_refseq_to_uniprot():
    """Given a list of RefSeq accessions, returns the mapping to UniProt."""

    def map_refseq_to_uniprot_helper():
        refseq_accs = get_all_proteins()
        ncbi_records = fetch_all_ncbi_protein_records()
        all_mappings = {}

        uniprot_mapping = refseq_to_uniprot(refseq_accs)
        print "Found %d mappings." % len(uniprot_mapping)
        for refseq_acc in tqdm(refseq_accs):
            all_mappings[refseq_acc] = list(uniprot_mapping.get(refseq_acc, {}))

        # For not mapped NP and YP accessions, get their WP accession and try
        # mapping again.
        not_mapped_refseqs = [refseq_acc for refseq_acc in refseq_accs
                              if not all_mappings[refseq_acc] and
                              (refseq_acc.startswith('NP') or
                               refseq_acc.startswith('YP'))]
        print ("Getting WP accs for %d NP/YP accs." % len(not_mapped_refseqs))

        wp_accession_mapping = {refseq:parse_wp_accession(ncbi_records[refseq])
                                 for refseq in not_mapped_refseqs}
        # Dump this for reference.
        wp_accession_mapping_json = os.path.join(DATA_DIR, 'np_yp_to_wp.json')
        json.dump(wp_accession_mapping, open(wp_accession_mapping_json, 'w'),
                  indent=4)

        uniprot_mapping = refseq_to_uniprot(wp_accession_mapping.values())
        print "Found %d additional mappings" % len(uniprot_mapping)
        for not_mapped_refseq in not_mapped_refseqs:
            wp = wp_accession_mapping[not_mapped_refseq]
            all_mappings[not_mapped_refseq] = list(uniprot_mapping.get(wp, {}))

        print "UniProt mapping size distribution:"
        c = Counter(map(len, all_mappings.values()))
        pprint(c)

        return all_mappings

    json_dump_file = os.path.join(DATA_DIR, "refseq_to_uniprot.json")
    if not os.path.isfile(json_dump_file):
        all_mappings = map_refseq_to_uniprot_helper()
        json.dump(all_mappings, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def batch_ncbi_taxonomy_id():
    """Writes NCBI taxonomy id for each RefSeq accession number to a file."""
    ncbi_records = fetch_all_ncbi_protein_records()
    tax_ids = {refseq:parse_ncbi_taxonomy_id(record)
               for refseq, record in ncbi_records.items()}
    json_dump_file = os.path.join(DATA_DIR, 'refseq_taxonomy_ids.json')
    json.dump(tax_ids, open(json_dump_file, 'w'), indent=4)

def batch_uniprot_taxonomy_id():
    """Writes NCBI taxonomy id for each UniProt accession number to a file."""
    mappings = map_refseq_to_uniprot()
    uniprot_records = fetch_all_uniprot_records()
    # Get all UniProt accessions and UniParc accessions.
    uniprot_accessions = [acc for accs in mappings.values() for acc in accs
                          if not acc.startswith('UPI')]
    tax_ids = {acc:parse_uniprot_tax_id(uniprot_records[acc])
               for acc in uniprot_accessions}
    json_dump_file = os.path.join(DATA_DIR, 'uniprot_taxonomy_ids.json')
    json.dump(tax_ids, open(json_dump_file, 'w'), indent=4)

def filter_same_taxon_accessions(mappings):
    """Filters mappings that have same NCBI taxonomy ID for both RefSeq and
    UniProt accession numbers.
    """
    ncbi_records = fetch_all_ncbi_protein_records()
    ambiguous_mappings = {k:v for k, v in mappings.items() if len(v) > 1}
    print "Number of ambiguous mappings:", len(ambiguous_mappings)
    same_taxon_mappings = {}

    for refseq_acc in tqdm(ambiguous_mappings):
        print refseq_acc, mappings[refseq_acc]
        refseq_tax_id = extract_ncbi_taxonomy_id(refseq_rec)
        uniprot_accs = [uniprot_acc for uniprot_acc in mappings[refseq_acc]
                        if get_uniprot_tax_id(uniprot_acc) == refseq_tax_id]
        if uniprot_acc:
            mappings[refseq_acc] = uniprot_accs
    return mappings

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
