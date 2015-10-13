"""Generates the file for linking CollecTF to UniProt.

The generated file is tab-separated and the format is as follows:

Uniprot Acc | Primary ID   | Opt 1 | Opt 2
Q9A724      | EXPREG_00001 | Collection of X YY binding sites | consensus seq.

where X=# of sites YY=TF name
"""

from base import bioutils
from base import models
from browse import dbxref

def generate_uniprot_dbxref():
    for TF_instance in models.TFInstance.objects.all():
        print TF_instance.uniprot_accession,
        print dbxref.to_uniprot_dbxref(TF_instance.TF_instance_id),
        site_ids = models.Curation_SiteInstance.objects.filter(
            site_type__in=['motif_associated', 'var_motif_associated'],
            curation__TF_instances=TF_instance).values_list(
                'site_instance', flat=True).distinct()
        sites = models.SiteInstance.objects.filter(site_id__in=site_ids)
        print "Collection of %d %s binding sites" % (sites.count(), TF_instance.TF),
        print bioutils.degenerate_consensus(
            bioutils.build_motif(site.seq for site in sites)),
            


def run():
    generate_uniprot_dbxref()
