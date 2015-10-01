"""
Given RefSeq accession numbers of the TF-instances in CollecTF, find their
UniProt identifiers.
"""

from collections import Counter
import csv
import json
import os
import pickle
import re
import xmltodict

from Bio import Entrez
from tqdm import tqdm
import uniprot

#from base import models

DATA_DIR = '/home/sefa/Dropbox/collectf/scripts/uniprot_migration/data'
#DATA_DIR = '/Users/sefa/Dropbox/collectf/scripts/uniprot_migration/data'
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

def parse_uniprot_record_review_status(record):
    """Returns true if the UniProt protein record is reviewed."""
    status = re.search(r'ID.*(Reviewed|Unreviewed);', record).group(1)
    assert status in ['Reviewed', 'Unreviewed']
    return status == 'Reviewed'

def does_uniprot_have_proteome(record):
    """Returns if the UniProt record is a match for a complete proteome."""
    return bool(re.search(r'DR   Proteomes;.*', record))

def is_refseq_accession(accession):
    """Returns true if the given accession number is from RefSeq."""
    prefixes = ['NP', 'YP', 'WP']
    return any(accession.startswith(prefix) for prefix in prefixes)

def is_np_or_yp_accession(accession):
    """Returns true if the accession is a NP or YP."""
    return accession.startswith('NP') or accession.startswith('YP')

def is_wp_accession(accession):
    """Returns true if the accession is a WP one."""
    assert is_refseq_accession(accession)
    return accession.startswith('WP')

def is_uniprot_accession(acc):
    """Returns true if the given accession is NOT RefSeq accession number."""
    return bool(re.match(
        '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',
        acc))

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
        uniprot_mapping = batch_refseq_to_uniprot()
        uniprot_accs = [acc for accs in uniprot_mapping.values()
                        for acc in accs if not is_uniparc_accession(acc)]

        all_records = {}
        print "Number of UniProt accession numbers:", len(uniprot_accs)

        for acc in tqdm(uniprot_accs):
            record, = fetch_uniprot_records([acc])
            all_records[acc] = record

        json.dump(all_records, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def refseq_to_uniprot(refseq_accession):
    """Maps the given RefSeq accession to UniProt ID."""
    return uniprot.map(refseq_accession, f='P_REFSEQ_AC', t='ACC')

def uniparc_to_uniprot(uniparc_accession):
    """Maps the given UniParc accession to set of UniProt accessions."""
    assert is_uniparc_accession(uniparc_accession)
    return uniprot.map(uniparc_accession, f='UPARC', t='ACC')

def batch_uniparc_to_uniprot(uniparc_accessions):
    """Maps all given UniParc accessions to UniProt accessions."""
    assert all(is_uniparc_accession(acc) for acc in uniparc_accessions)
    uniprot_map = uniprot.map(uniparc_accessions, f='UPARC', t='ACC')
    # Dump to file for reference.
    json_file = os.path.join(DATA_DIR, 'uniparc_to_uniprot.json')
    json.dump({k: list(v) for k, v in uniprot_map.items()},
              open(json_file, 'w'), indent=4)
    return uniprot_map

def batch_refseq_to_uniprot():
    """Given a list of RefSeq accessions, returns the mapping to UniProt."""
    json_dump_file = os.path.join(DATA_DIR, "simple_refseq_to_uniprot.json")
    if not os.path.isfile(json_dump_file):
        all_mapping = {}
        refseq_accs = get_all_proteins()
        uniprot_mapping = refseq_to_uniprot(refseq_accs)
        for refseq_acc in tqdm(refseq_accs):
            all_mapping[refseq_acc] = list(uniprot_mapping.get(refseq_acc, {}))

        # Some NP/YP accessions are not matched. Map them from their WP
        # accessions.
        map_not_matched_np_and_yps(all_mapping)
        # UniProt API maps some accessions to UniParc. Get their UniProt
        # accessions.
        map_uniparcs_back_to_uniprot(all_mapping)

        json.dump(all_mapping, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def map_not_matched_np_and_yps(mapping):
    """Maps not-matched NP/YP accessions to WP and attempts UniProt mapping."""
    not_mapped_refseqs = [refseq_acc for refseq_acc in mapping
                          if not mapping[refseq_acc] and
                          is_np_or_yp_accession(refseq_acc)]
    print "Getting WP accs for %d NP/YP accs." % len(not_mapped_refseqs)
    ncbi_records = fetch_all_ncbi_protein_records()
    wp_accession_mapping = {refseq:parse_wp_accession(ncbi_records[refseq])
                            for refseq in not_mapped_refseqs}
    uniprot_mapping = refseq_to_uniprot(wp_accession_mapping.values())
    print "Found %d additional matchings" % len(uniprot_mapping)
    for not_mapped_refseq in not_mapped_refseqs:
        wp = wp_accession_mapping[not_mapped_refseq]
        mapping[not_mapped_refseq] = list(uniprot_mapping.get(wp, {}))

    # Dump NP/YP -> WP mapping for reference.
    wp_accession_mapping_json = os.path.join(DATA_DIR, 'np_yp_to_wp.json')
    json.dump(wp_accession_mapping, open(wp_accession_mapping_json, 'w'),
              indent=4)

def map_uniparcs_back_to_uniprot(mapping):
    """Maps UniParc accessions back to set of UniProt accessions.

    UniProt maps some RefSeq accessions to UniParc accessions although UniProt
    is selected for target DB, not UniParc. For such cases, the solution is to
    map UniParc accession back to UniProt accessions.
    """
    # Collect all UniParc accessions for batch UniProt mapping.
    uniparc_accs = [acc for accs in mapping.values() for acc in accs
                    if is_uniparc_accession(acc)]
    uniparc_to_uniprot_map = batch_uniparc_to_uniprot(uniparc_accs)
    for refseq in tqdm(mapping):
        for acc in mapping[refseq]:
            # Replace UniParc accessions with corresponding UniProt ones.
            if is_uniparc_accession(acc):
                mapping[refseq].remove(acc)
                assert acc not in mapping[refseq]
                mapping[refseq].extend(uniparc_to_uniprot_map.get(acc, {}))

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
        mapping = batch_refseq_to_uniprot()
        uniprot_records = fetch_all_uniprot_records()
        # Get all UniProt accessions and UniParc accessions.
        uniprot_accessions = [acc for accs in mapping.values() for acc in accs
                              if not acc.startswith('UPI')]
        tax_ids = {acc: parse_uniprot_tax_id(uniprot_records[acc])
                   for acc in uniprot_accessions}
        json.dump(tax_ids, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def map_to_acc_with_same_taxon(mapping, mapping_notes):
    """Filters mapping that have same NCBI taxonomy ID for both RefSeq and
    UniProt accession numbers.
    """
    ncbi_records = fetch_all_ncbi_protein_records()
    ncbi_taxonomy_ids = batch_ncbi_taxonomy_id()
    uniprot_taxonomy_ids = batch_uniprot_taxonomy_id()
    same_taxon_mapping = {}

    for refseq_acc in tqdm(mapping):
        if len(mapping[refseq_acc]) <= 1:
            continue
        refseq_tax_id = ncbi_taxonomy_ids[refseq_acc]
        uniprot_accs = [acc for acc in mapping[refseq_acc]
                        if uniprot_taxonomy_ids[acc] == refseq_tax_id]
        if uniprot_accs:
            mapping[refseq_acc] = [uniprot_accs[0]]
            mapping_notes[refseq_acc] = "taxonomy id match."

def batch_uniprot_review_status():
    """Checks UniProt records for their review status."""
    uniprot_records = fetch_all_uniprot_records()
    status = {acc: parse_uniprot_record_review_status(rec)
              for acc, rec in uniprot_records.items()}
    json_dump_file = os.path.join(DATA_DIR, 'uniprot_review_status.json')
    json.dump(status, open(json_dump_file, 'w'), indent=4)
    return status

def map_to_reviewed_accessions(mapping, mapping_notes):
    """Resolves ambiguous matchings with reviewed UniProt records."""
    uniprot_review_status = batch_uniprot_review_status()
    for refseq, uniprots in mapping.items():
        if len(uniprots) <= 1:
            continue
        reviewed_uniprots = [uniprot for uniprot in uniprots
                             if uniprot_review_status[uniprot]]
        if reviewed_uniprots:
            mapping[refseq] = [reviewed_uniprots[0]]
            mapping_notes[refseq] = "reviewed UniProt record"

def parse_proteome_file():
    """Parses proteome file.

    For each proteome, parses
    - its taxonomy id (taxid)
    - if it is representative for the proteome cluster
    - if it belongs to reference species/strain.
    """
    json_dump_file = os.path.join(DATA_DIR, 'proteomes.json')
    if not os.path.isfile(json_dump_file):
        proteome_file = os.path.join(DATA_DIR, 'proteomes_2015_09.xml')
        with open(proteome_file) as f:
            obj = xmltodict.parse(f.read())
        proteomes = []
        for proteome in obj['proteomes']['proteome']:
            proteomes.append(
                {'taxid': proteome['taxonomy'],
                 'is_reference': proteome['is_reference_proteome'] == 'true',
                 'is_representative':
                     proteome['is_representative_proteome'] == 'true'})
        json.dump(proteomes, open(json_dump_file, 'w'), indent=4)

    return json.load(open(json_dump_file))

def map_to_accessions_with_proteome(mapping, mapping_notes):
    """Resolves mapping using proteome data.

    For each candidate UniProt accession, pick
    - if its taxonomy id is a match for reference/representative proteome or
    - if its record has match for a proteome record.
    """
    # Check if taxonomy id is a match for a reference/representative proteome
    proteomes = parse_proteome_file()
    uniprot_taxonomy_ids = batch_uniprot_taxonomy_id()
    for refseq, uniprots in mapping.items():
        if len(uniprots) <= 1:
            continue
        for uniprot in uniprots:
            uniprot_taxid = uniprot_taxonomy_ids[uniprot]
            if any(proteome['taxid'] == uniprot_taxid and
                   (proteome['is_reference'] or proteome['is_representative'])
                   for proteome in proteomes):
                mapping[refseq] = [uniprot]
                mapping_notes[refseq] = (
                    "taxon match to reference/representative proteome")
                break
    mapping_stats(mapping)

    # Check if the UniProt record has match for a proteome record.
    uniprot_records = fetch_all_uniprot_records()
    for refseq, uniprots in mapping.items():
        if len(uniprots) <= 1:
            continue
        for uniprot in uniprots:
            if does_uniprot_have_proteome(uniprot_records[uniprot]):
                mapping[refseq] = [uniprot]
                mapping_notes[refseq] = "UniProt record matches to a proteome"
                break

def mock_get_all_proteins():
    proteins =[u'NP_799324', u'NP_231738', u'YP_006516164',
               u'YP_006516595', u'NP_417816']
    return proteins

def get_ambiguous_matchings(mapping):
    """Returns matchings with more than one UniProt accessions."""
    return [accs for accs in mapping.values() if len(accs) > 1]

def get_missing_matchings(mapping):
    """Retuns matchings with no UniProt accessions."""
    return [accs for accs in mapping.values() if not accs]

def mapping_stats(mapping):
    print "num missing matchings:", len(get_missing_matchings(mapping))
    print "num ambigous matchings:", len(get_ambiguous_matchings(mapping))
    print "//"

def resolve_mapping():
    """Gets initial mapping and resolves  to multiple UniProt accessions

    - Checks if any of the identified UniProt accessions have the same taxonomy
      id as the RefSeq accession.
    - Checks if any of the UniProt accessions have 'reviewed' status and uses it
      if any.
    """
    # Get all mappings
    mapping = batch_refseq_to_uniprot()
    mapping_stats(mapping)

    mapping_notes = {}
    for refseq, matches in mapping.items():
        if len(matches) == 1:
            mapping_notes[refseq] = \
                "single UniProt accession, no resolution required"
    # Resolve ambiguous matching using taxonomy. For each RefSeq accession, if
    # there are more than one matching UniProt accession, pick the one that has
    # the same taxonomy ID as RefSeq accession has.
    map_to_acc_with_same_taxon(mapping, mapping_notes)
    mapping_stats(mapping)

    # For each RefSeq accession with more than one matching UniProt accessions,
    # pick the first reviewed UniProt accession.
    map_to_reviewed_accessions(mapping, mapping_notes)
    mapping_stats(mapping)

    # Resolve using proteome data
    map_to_accessions_with_proteome(mapping, mapping_notes)
    mapping_stats(mapping)

    json_dump_file = os.path.join(DATA_DIR, 'resolved_refseq_to_uniprot.json')
    json.dump(mapping, open(json_dump_file, 'w'), indent=4)

    csv_file = os.path.join(DATA_DIR, 'resolved_refseq_to_uniprot.csv')
    with open(csv_file, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['RefSeq accession', 'UniProt accession', 'Notes'])
        for refseq in mapping:
            if mapping[refseq]:
                assert len(mapping[refseq]) == 1
                writer.writerow([refseq, mapping[refseq][0],
                                 mapping_notes[refseq]])
            else:
                writer.writerow([refseq, "N/A", "no matching found"])
