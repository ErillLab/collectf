"""The module for generating CollecTF release file."""

from base import models
from browse import motif_report
from browse import export

def make_new_release():
    """Generates the release file containing all binding sites in CollecTF."""
    all_csis = models.Curation_SiteInstance.objects.all()
    reports = motif_report.make_distinct_reports(all_csis)
    # Get the list of meta sites where each meta site is a list of
    # Curation_SiteInstance ids
    meta_site_id_lists = [[csi.id for csi in ms.cur_site_insts]
                          for report in reports
                          for ms in report.get_meta_sites()]

    meta_sites = [models.Curation_SiteInstance.objects.filter(pk__in=ms_id_list)
                  for ms_id_list in meta_site_id_lists]

    tsv_raw = export.export_tsv(meta_sites)
    with open('collectf_export.tsv', 'w') as f:
        f.write(tsv_raw)

def run():
    """Entry point for the script."""
    make_new_release()
