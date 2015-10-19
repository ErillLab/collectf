"""This module generates the file containing motifs in CollecTF."""

import os
from tqdm import tqdm

from base import models
from base.bioutils import build_motif
from browse.static_reports import get_static_reports
from collectf import settings

def generate_motif_file():
    """Generates motif file containing all motifs for each TF-instance."""
    export_file = os.path.join(settings.STATICFILES_DIRS[0], 'collectf_meme.cm')
    with open(export_file, 'w') as f:
        for TF_instance in tqdm(models.TFInstance.objects.all()):
            reports, _ = get_static_reports(
                'TF_instance_%s' % TF_instance.TF_instance_id)
            if reports:
                motif_sites = reports[0].get_single_motif()
                if len(motif_sites) < 8:
                    continue
                motif = build_motif(motif_sites)
                f.write('>%s\n' % TF_instance.uniprot_accession)
                f.write('A| %s\n' % ' '.join(map(str, motif.counts['A'])))
                f.write('C| %s\n' % ' '.join(map(str, motif.counts['C'])))
                f.write('G| %s\n' % ' '.join(map(str, motif.counts['G'])))
                f.write('T| %s\n' % ' '.join(map(str, motif.counts['T'])))

def run():
    """Entry point for the script."""
    generate_motif_file()
