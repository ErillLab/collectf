from collectfapp import models
from collectfapp import bioutils
from browseapp import search
from browseapp import view_results
from browseapp import export
from collectf import settings
import os

def make_new_release():
    """Make a release file of the entire database."""
    all_csis = models.Curation_SiteInstance.objects.all()
    motif_csi_ids = all_csis.filter(is_motif_associated=True).values_list('pk', flat=True)
    non_motif_csi_ids = all_csis.filter(is_motif_associated=False).values_list('pk', flat=True)
    motif_csi_ids = list(motif_csi_ids)
    non_motif_csi_ids = list(non_motif_csi_ids)
    template = view_results.prepare_results(motif_csi_ids, non_motif_csi_ids, True)
    reports = template['reports'] # each report is all binding sites for (TF, species) pair
    # each ms is a list of curation site instance objects
    meta_sites = [ms for report in reports for ms in report['meta_sites'].values()]
    # convert the list of meta sites into queryset
    meta_sites = [models.Curation_SiteInstance.objects.filter(pk__in=[csi.pk for csi in ms])
                  for ms in meta_sites]
    tsv_raw = export.export_tsv_raw(meta_sites)
    with open(os.path.join(settings.MEDIA_ROOT, 'collectf_export.tsv'), 'w') as f:
        f.write(tsv_raw)


    
def run():
    make_new_release()
