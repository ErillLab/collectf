"""Generate release file"""

import os
from base import models
from browse import motif_report
from browse import export

def make_new_release():
    """Make a release file of the entire database.
    """
    all_csis = models.Curation_SiteInstance.objects.all()
    motif_assoc = all_csis.filter(site_type='motif_associated')
    non_moitf_assoc = all_csis.filter(site_type='non_motif_associated')
    # group sites by TF and genome
    vals = motif_assoc.values_list('curation__TF_instances',
                                   'site_instance__genome')\
                      .distinct()\
                      .order_by('curation__TF__name',
                                'site_instance__genome')

    meta_sites = []
    for (TF_instances, genome) in vals:
        meta_sites.append(models.Curation_SiteInstance.objects.filter(
            curation__TF_instances=TF_instances,
            site_instance__genome=genome))

    tsv_raw = export.export_tsv(meta_sites)
    with open('collectf_export.tsv', 'w') as f:
        f.write(tsv_raw)



def run():
    make_new_release()
