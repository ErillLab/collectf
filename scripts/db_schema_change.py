"""
Scripts for CollecTF database schema changes. Intented to run once.
"""

import csv
import pickle
import re
from tqdm import tqdm

from Bio import Entrez

from base import models

Entrez.email = 'sefa1@umbc.edu'

def move_TF_from_curation_to_TF_instance_table():
    """Moves relationship between curation and TF to TF and TF-instance."""
    for curation in tqdm(models.Curation.objects.all()):
        for TF_instance in curation.TF_instances.all():
            TF_instance.TF = curation.TF
            TF_instance.save()

def fetch_ncbi_protein_record(accession):
    """Fetches the protein record from NCBI, for the given accession number."""
    handle = Entrez.efetch(db='protein', id=accession, retmode='xml')
    records = Entrez.read(handle)
    return records[0]

def parse_wp_accession(rec):
    """Returns a WP accession for a given NP/YP accession."""
    contig = rec['GBSeq_contig']
    wp_acc = re.match(r'join\((WP_\d+.\d):\d..\d+\)', contig).group(1)
    return wp_acc.split('.')[0]

def refseq_to_wp(acc):
    """Maps given RefSeq (NP or YP) accession number to WP one."""
    no_move_list = ['NP_059642'] # Don't move these accessions to WP.
    if not (acc.startswith('WP') or acc in no_move_list):
        rec = fetch_ncbi_protein_record(acc)
        wp = parse_wp_accession(rec)
        print acc, wp
        return wp
    return acc

def add_refseq_accession():
    """Populates WP accession numbers on 'refseq_accession' field."""
    for TF_instance in models.TFInstance.objects.all():
        print TF_instance.protein_accession
        TF_instance.refseq_accession = refseq_to_wp(TF_instance.protein_accession)
        if TF_instance.refseq_accession != TF_instance.protein_accession:
            TF_instance.notes = (
                "Original RefSeq accession: %s" % TF_instance.protein_accession)
        TF_instance.save()

def merge_same_wps():
    """Merges TF instances with the same WP accession number."""
    wp_accessions = models.TFInstance.objects.values_list(
        'refseq_accession', flat=True).distinct()
    for wp_acc in wp_accessions:
        tf_instances = models.TFInstance.objects.filter(refseq_accession=wp_acc)
        if tf_instances.count() == 1:
            continue
        # Get the first one and merge the rest
        primary_tf_instance = tf_instances[0]
        print primary_tf_instance.protein_accession, '--',
        for obsolete_tf_instance in tf_instances[1:]:
            print obsolete_tf_instance.protein_accession,
            for curation in models.Curation.objects.filter(
                    TF_instances=obsolete_tf_instance):
                print curation.curation_id,
                curation.TF_instances.remove(obsolete_tf_instance)
                curation.TF_instances.add(primary_tf_instance)
            assert not models.Curation.objects.filter(
                TF_instances=obsolete_tf_instance)
        
            # Add note about the merge on primary TF-instance record
            primary_tf_instance.notes += (
                '\n\n%s merged to this record.\n' %
                obsolete_tf_instance.protein_accession)
            primary_tf_instance.save()
            obsolete_tf_instance.delete()
        print ''

def add_uniprot_accessions():
    """Populates UniProt field on TF-instance table."""
    csv_file = 'scripts/uniprot_migration/data/resolved_refseq_to_uniprot.csv'
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        mapping = {row['RefSeq accession']: row['UniProt accession']
                   for row in reader}
    for TF_instance in models.TFInstance.objects.all():
        print TF_instance.protein_accession, mapping[TF_instance.protein_accession]
        TF_instance.uniprot_accession = mapping[TF_instance.protein_accession]
        TF_instance.save()

def dump_curation_TF_instance():
    """Dumps Curation-TF_instance table to a a file."""
    rows = []
    for curation in models.Curation.objects.all():
        for uniprot, refseq, name in curation.TF_instances.values_list(
                'uniprot_accession', 'refseq_accession', 'name'):
            print curation.curation_id, uniprot, refseq, name
            rows.append((curation.curation_id, uniprot))
        print '---'

    pickle.dump(rows, open('curation_TF_instance_dump.pkl', 'w'))

def load_curation_TF_instance():
    """Load Curation-TF_instance M2M table from given pickle file."""
    rows = pickle.load(open('curation_TF_instance_dump.pkl'))
    for row in rows:
        curation_id, uniprot = row
        curation = models.Curation.objects.get(curation_id=curation_id)
        TF_instance = models.TFInstance.objects.get(uniprot_accession=uniprot)
        curation.TF_instances.add(TF_instance)
        curation.save()

def dump_TF_instance():
    """Dumps TF-instance table to a pickle file."""
    rows = []
    for tfi in models.TFInstance.objects.all():
        print tfi.protein_accession,
        print tfi.refseq_accession,
        print tfi.uniprot_accession,
        print tfi.description,
        print tfi.TF.name,
        print tfi.notes
        rows.append((tfi.protein_accession,
                     tfi.refseq_accession,
                     tfi.uniprot_accession,
                     tfi.description,
                     tfi.TF.name,
                     tfi.notes))

    pickle.dump(rows, open('TF_instance_dump.pkl', 'w'))

def load_TF_instance():
    """Loads TF instance table from given pickle file."""
    rows = pickle.load(open('TF_instance_dump.pkl'))
    for row in rows:
        TF = models.TF.objects.get(name=row[4])
        TF_instance = models.TFInstance(
            refseq_accession=row[1],
            uniprot_accession=row[2],
            description=row[3],
            TF=TF,
            notes=row[5])
        TF_instance.save()


# TODO(sefa): update TF_instance descriptions
def udpate_TF_instance_description():
    pass

def run():
    """Entry point for the script."""
    #move_TF_from_curation_to_TF_instance_table()
    #refseq_to_uniprot_accession()
    #move_to_wp()
    #add_refseq_accessions()
    #merge_same_wps()
    #add_uniprot_accessions()
    #dump_curation_TF_instance()
    #dump_TF_instance()
    #load_TF_instance()
    load_curation_TF_instance()
