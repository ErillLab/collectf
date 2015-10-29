"""Generates the file for linking CollecTF to UniProt.

The generated file is tab-separated and the format is as follows:

Uniprot Acc | Primary ID   | Opt 1 | Opt 2
Q9A724      | EXPREG_00001 | Collection of X YY binding sites | consensus seq.

where X=# of sites YY=TF name
"""

import os

from tqdm import tqdm

from base import models
from base.bioutils import build_motif
from base.bioutils import degenerate_consensus
from browse import dbxref
from browse import motif_report
from collectf import settings

def generate_uniprot_dbxref():
    export_file = os.path.join(
        settings.STATIC_ROOT, 'uniprot_dbxref.txt')
    with open(export_file, 'w') as f:
        for TF_instance in models.TFInstance.objects.all():
            print TF_instance
            curation_site_instances = models.Curation_SiteInstance.objects.filter(
                curation__TF_instances=TF_instance)
            reports = motif_report.make_reports(curation_site_instances)
            if reports:
                f.write('\t'.join(
                    [TF_instance.uniprot_accession,
                     dbxref.to_uniprot_dbxref(TF_instance.TF_instance_id)]))
                f.write('\n')
def run():
    generate_uniprot_dbxref()
