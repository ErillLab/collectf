"""Generate release file"""

import os
from base import models
from browse import motif_report
from browse import export

def make_new_release():
    """Make a release file of the entire database.
    """
    all_csis = models.Curation_SiteInstance.objects.all()
    reports = motif_report.make_reports(all_csis)
    # get delegate sites for each meta site
    meta_sites_query_set = []
    for report in reports:
        for meta_site in report.get_meta_sites():
            qset = models.Curation_SiteInstance.objects.filter(
                pk__in=[csi.pk for csi in meta_site.cur_site_insts])
            meta_sites_query_set.append(qset)


    tsv_raw = export.export_tsv(meta_sites_query_set)
    with open('collectf_export.tsv', 'w') as f:
        f.write(tsv_raw)



def run():
    make_new_release()
