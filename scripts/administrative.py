"""Collection of code snippets for internal use"""

from collections import defaultdict
import csv
import pandas as pd
import time

from tqdm import tqdm

import base
from curate import create_object
from base import models

def get_TFs():
    """Gets the list of TFs and their descriptions."""
    tfs = models.TF.objects.all()
    for tf in tfs:
        print tf.name, "\t", tf.description

def add_pubs_from_csv(filename):
    """Batch publication submission"""
    with open(filename) as f:
        df = [line.strip().split(',') for line in f.readlines()[1:]]

    for row in df:
        print row[0], row[1], row[2]
        try:
            p = models.Publication.objects.get(pmid=row[0])
        except models.Publication.DoesNotExist:
            pubrec = base.bioutils.get_pubmed(row[0])
            cd = dict(pmid=row[0],
                      reported_TF=row[1].strip(),
                      reported_species=row[2].strip(),
                      contains_promoter_data=False,
                      contains_expression_data=False,
                      submission_notes="")

            time.sleep(2)
            p = create_object.make_pub(pubrec, cd)
            p.save()
            print p
    return df

def validate_curations():
    sefa = models.Curator.objects.get(user__username='sefa1')
    for curation in models.Curation.objects.all():
        if curation.curator.user.username in ['dinara', 'ivanerill', 'ewhite5@umbc.edu']:
            curation.validated_by = sefa
            curation.save()
            print "Curation", curation.pk, "validated."

        elif int(curation.pk) < 554:
            curation.validated_by = sefa
            curation.save()
            print "Curation", curation.pk, "validated."

def pub_analysis():
    pubs = models.Publication.objects.all()
    s = pd.Series([pub.publication_date.split()[0] for pub in pubs], name='date')
    s.to_csv("pub_dates.csv", index=False, header=True)

def site_analysis():
    sites = models.Curation_SiteInstance.objects.filter(site_type='motif_associated')
    s = pd.Series([site.curation.created or "2013-08-06 16:24:09+00:00"
                   for site in sites], name='time')
    s.to_csv("site_dates.csv", index=False, header=True)

def list_pubs():
    pubs = models.Publication.objects.all()
    with open("pub_list.csv", 'w') as f:
        f.write('\t'.join(['PMID', 'Journal', 'Completed']) + '\n')
        for pub in pubs:
            f.write('%s\t%s\t%s\n' % (pub.pmid, pub.journal, pub.curation_complete))

def collectf_tf_instance_summary():
    """Lists all TF instances and their TFs and the curation IDs into a file."""
    with open('scripts/data/tf_instance_summary.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['TF instance', 'TF', 'Curations using the TF instance'])
        tf_instances = models.TFInstance.objects.all()
        for tf_instance in tqdm(tf_instances):
            # Get curations using this TF instance.
            curations = models.Curation.objects.filter(TF_instances=tf_instance)
            tfs = set(curation.TF for curation in curations)
            writer.writerow([tf_instance.protein_accession,
                             ', '.join([tf.name for tf in tfs]),
                             ', '.join([str(curation.curation_id)
                                        for curation in curations])])

def run():
    collectf_tf_instance_summary()
