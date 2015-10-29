"""This module generates the file containing motifs in CollecTF."""

import os

from base import models
from base.bioutils import build_motif
from browse import dbxref
from browse import motif_report
from collectf import settings

def generate_motif_file():
    """Generates motif file containing all motifs for each TF-instance."""
    export_file = os.path.join(settings.STATIC_ROOT, 'collectf_meme.cm')
    with open(export_file, 'w') as f:
        for TF_instance in models.TFInstance.objects.all():
            print TF_instance
            curation_site_instances = models.Curation_SiteInstance.objects.filter(
                curation__TF_instances=TF_instance)
            reports = motif_report.make_reports(curation_site_instances)
            if reports:
                sites = reports[0].align_sites()
                if len(sites) < 10:
                    continue
                motif = build_motif(sites)
                f.write('>%s %s_%s\n' % (
                    dbxref.to_uniprot_dbxref(TF_instance.TF_instance_id),
                    reports[0].TF_name,
                    reports[0].species_name.replace(' ', '_')))
                f.write('%s\n' % '\t'.join(map(str, motif.counts['A'])))
                f.write('%s\n' % '\t'.join(map(str, motif.counts['C'])))
                f.write('%s\n' % '\t'.join(map(str, motif.counts['G'])))
                f.write('%s\n' % '\t'.join(map(str, motif.counts['T'])))

def run():
    """Entry point for the script."""
    generate_motif_file()
