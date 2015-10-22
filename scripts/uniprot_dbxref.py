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
from browse.static_reports import get_static_reports
from collectf import settings

def generate_uniprot_dbxref():
    export_file = os.path.join(
        settings.STATICFILES_DIRS[0], 'uniprot_dbxref.txt')
    with open(export_file, 'w') as f:
        for TF_instance in tqdm(models.TFInstance.objects.all()):
            reports, _ = get_static_reports(
                'TF_instance_%s' % TF_instance.TF_instance_id)
            if reports:
                motif_sites = reports[0].get_single_motif()
                f.write('\t'.join(
                    [TF_instance.uniprot_accession,
                     dbxref.to_uniprot_dbxref(TF_instance.TF_instance_id),
                     "Collection of %d %s binding sites" %
                     (len(motif_sites), TF_instance.TF.name),
                     degenerate_consensus(build_motif(motif_sites))]))
                f.write('\n')
def run():
    generate_uniprot_dbxref()
