"""The module for generating CollecTF release file."""

import os
from zipfile import ZipFile

from django.conf import settings

from core import models
from browse import motif_report
from browse import export


def db_dump():
    """Generates the release file containing all binding sites in CollecTF."""
    curation_site_instances = models.Curation_SiteInstance.objects.all()
    reports = motif_report.build_motif_reports(curation_site_instances)
    meta_sites = [meta_site
                  for report in reports
                  for meta_site in report.meta_sites]
    tsv_raw = export.export_tsv(meta_sites)
    export_file = os.path.join(
        settings.STATICFILES_DIRS[0], 'collectf_export.tsv')
    with open(export_file, 'w') as f:
        f.write(tsv_raw)


def all_alignments():
    """Generates aligned binding moitfs."""
    curation_site_instances = models.Curation_SiteInstance.objects.all()
    reports = motif_report.build_motif_reports(curation_site_instances)
    zip_file_name = os.path.join(settings.STATICFILES_DIRS[0],
                                 'collectf_alignments.zip')
    with ZipFile(zip_file_name, 'w') as zip:
        for report in reports:
            if len(report.aligned_sites) >= 8:
                tf = report.TF_name
                species = report.short_species_name.replace('.', '_')
                print tf, species
                filename = '%s_%s.txt' % (tf, species)
                zip.writestr(filename, '\n'.join(report.aligned_sites))


def run():
    """Entry point for the script."""
    db_dump()
