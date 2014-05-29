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
    import pandas as pd
    df = pd.read_csv(filename, sep=',') # read as Pandas data frame
    for i, row in df.iterrows():
        print row.PMID, row.TF, row.Strain
        try:
            p = models.Publication.objects.get(pmid=row.PMID)
        except models.Publication.DoesNotExist:
            pubrec = base.bioutils.get_pubmed(row.PMID)
            cd = dict(pmid=row.PMID,
                      reported_TF=row.TF.strip(),
                      reported_species=row.Strain.strip(),
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
