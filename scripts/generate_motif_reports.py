# This script generates motif reports for all TFs, species and experimental
# technique levels. The generated motif reports are saved and used for browse
# pages in CollecTF without accessing database and computing motif reports for
# the query TFs, species or experimental techniques.

from tqdm import tqdm

from base import models
from browse import motif_report

def motif_reports_by_tf_family():
    """Generates motif reports for each TF-family."""
    TF_families = models.TFFamily.objects.all()
    for TF_family in tqdm(TF_families):
        TFs = models.TF.objects.filter(family=TF_family)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF__in=TFs)
        reports = motif_report.make_reports(curation_site_instances)

def run():
    motif_reports_by_tf_family()


