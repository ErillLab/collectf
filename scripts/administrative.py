"""Collection of code snippets for internal use"""

import base
from curate import create_object
from base import models
import time

def get_TFs():
    """Get the list of TFs and their descriptions"""
    TFs = models.TF.objects.all()
    for TF in TFs:
        print TF.name, "\t", TF.description

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

def run():
    #get_TFs()
    add_pubs_from_csv("/home/sefa/Desktop/Book1.csv")
