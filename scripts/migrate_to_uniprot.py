"""
Given RefSeq accession numbers of the TF-instances in CollecTF, find their
UniProt identifiers.
"""

from collections import Counter
import json
import os
import pickle
import re
import xmltodict

from Bio import Entrez
from tqdm import tqdm
import uniprot

#from base import models

DATA_DIR = '/home/sefa/Dropbox/collectf/scripts/data'
#DATA_DIR = '/Users/sefa/Dropbox/collectf/scripts/data'
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

def uniprot_record_review_status(record):
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

        for acc in tqdm(uniprot_accs):
            record, = fetch_uniprot_records([acc])
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
        print "Getting WP accs for %d NP/YP accs." % len(not_mapped_refseqs)

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

    json_dump_file = os.path.join(DATA_DIR, "simple_refseq_to_uniprot.json")
    if not os.path.isfile(json_dump_file):
        all_mappings = map_refseq_to_uniprot_helper()
        json.dump(all_mappings, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def batch_ncbi_taxonomy_id():
    """Writes NCBI taxonomy id for each RefSeq accession number to a file."""
    json_dump_file = os.path.join(DATA_DIR, 'refseq_taxonomy_ids.json')
    if not os.path.isfile(json_dump_file):
        ncbi_records = fetch_all_ncbi_protein_records()
        tax_ids = {refseq:parse_ncbi_taxonomy_id(record)
                   for refseq, record in ncbi_records.items()}
        json.dump(tax_ids, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def batch_uniprot_taxonomy_id():
    """Writes NCBI taxonomy id for each UniProt accession number to a file."""
    json_dump_file = os.path.join(DATA_DIR, 'uniprot_taxonomy_ids.json')
    if not os.path.isfile(json_dump_file):
        mappings = map_refseq_to_uniprot()
        uniprot_records = fetch_all_uniprot_records()
        # Get all UniProt accessions and UniParc accessions.
        uniprot_accessions = [acc for accs in mappings.values() for acc in accs
                              if not acc.startswith('UPI')]
        tax_ids = {acc: parse_uniprot_tax_id(uniprot_records[acc])
                   for acc in uniprot_accessions}
        json.dump(tax_ids, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def same_taxon_accessions():
    """Filters mappings that have same NCBI taxonomy ID for both RefSeq and
    UniProt accession numbers.
    """
    json_dump_file = os.path.join(DATA_DIR, 'same_taxon_accessions.json')
    if not os.path.isfile(json_dump_file):
        mappings = map_refseq_to_uniprot()
        ncbi_records = fetch_all_ncbi_protein_records()
        ncbi_taxonomy_ids = batch_ncbi_taxonomy_id()
        uniprot_taxonomy_ids = batch_uniprot_taxonomy_id()
        same_taxon_mappings = {}

        for refseq_acc in tqdm(mappings):
            refseq_tax_id = ncbi_taxonomy_ids[refseq_acc]
            uniprot_accs = [acc for acc in mappings[refseq_acc]
                            if not is_uniparc_accession(acc) and
                            uniprot_taxonomy_ids[acc] == refseq_tax_id]
            if uniprot_accs:
                same_taxon_mappings[refseq_acc] = uniprot_accs

        json.dump(same_taxon_mappings, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def batch_uniprot_review_status():
    """Checks UniProt records for their review status."""
    uniprot_records = fetch_all_uniprot_records()
    status = {acc: uniprot_record_review_status(rec)
              for acc, rec in uniprot_records.items()}
    json_dump_file = os.path.join(DATA_DIR, 'uniprot_review_status.json')
    json.dump(status, open(json_dump_file, 'w'), indent=4)
    return status

def parse_proteome_file():
    """Parses proteome file.

    For each proteome, parses
    - its taxonomy id (taxid)
    - if it is representative for the proteome cluster
    - if it belongs to reference species/strain.
    """
    proteome_file = os.path.join(DATA_DIR, 'proteomes_2015_09.xml')
    with open(proteome_file) as f:
        obj = xmltodict.parse(f.read())
    proteomes = []
    for proteome in obj['proteomes']['proteome']:
        proteomes.append(
            {'taxid': proteome['taxonomy'],
             'is_reference': proteome['is_reference_proteome'],
             'is_representative': proteome['is_representative_proteome']})
    json_dump_file = os.path.join(DATA_DIR, 'proteomes.json')
    json.dump(proteomes, open(json_dump_file, 'w'), indent=4)

def mock_get_all_proteins():
    proteins =[u'NP_799324', u'NP_231738', u'YP_006516164',
               u'YP_006516595', u'NP_417816']
    return proteins

def migrate_to_uniprot():
    proteins = get_all_proteins()
    all_mappings = refseq_to_uniprot_batch(proteins)
    pickle.dump(all_mappings, open('all_mappings.pkl', 'w'))

def get_ambiguous_mappings(mappings):

    """Returns mappings with more than one UniProt accessions."""
    return [accs for accs in mappings.values() if len(accs) > 1]

def get_missing_mappings(mappings):
    """Retuns mappings with no UniProt accessions."""
    return [accs for accs in mappings.values() if not accs]

def mappings_stats(mappings):
    print "num missing mappings:", len(get_missing_mappings(mappings))
    print "num ambigous mappings:", len(get_ambiguous_mappings(mappings))
    print "//"

def resolve_mappings():
    """Gets initial mappings and resolves ambiguous mappings.

    - Checks if any of the identified UniProt accessions have the same taxonomy
      id as the RefSeq accession.
    - Checks if any of the UniProt accessions have 'reviewed' status and uses it
      if any.
    """

    # Find all mappings from RefSeq to UniProt.
    mappings = map_refseq_to_uniprot()
    mappings_stats(mappings)

    # Resolve ambiguous mappings using taxonomy.
    same_taxon_mappings = same_taxon_accessions()
    for refseq, uniprots in mappings.items():
        if len(uniprots) > 1 and len(same_taxon_mappings.get(refseq, [])) == 1:
            mappings[refseq] = same_taxon_mappings[refseq]
    mappings_stats(mappings)

    # Resolve ambiguous mappings with reviewed UniProt records.
    uniprot_review_status = batch_uniprot_review_status()
    for refseq, uniprots in mappings.items():
        if len(uniprots) > 1:
            reviewed_uniprots = [uniprot for uniprot in uniprots
                                 if uniprot_review_status[uniprot]]
            if len(reviewed_uniprots) > 0:
                mappings[refseq] = [reviewed_uniprots[0]]
    mappings_stats(mappings)


    json_dump_file = os.path.join(DATA_DIR, 'resolved_refseq_to_uniprot.json')
    json.dump(mappings, open(json_dump_file, 'w'), indent=4)
