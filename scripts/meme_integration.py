"""This module generates the file containing motifs in CollecTF."""

from tqdm import tqdm

from base import models
from browse.export import export_PSFM

def generate_motif_file():
    """Generates motif file containing all motifs for each TF-instance."""
    for TF_instance in tqdm(models.TFInstance.objects.all()):
        reports, _ = get_static_reports(
            'tf_instance_%s' % TF_instance.TF_instance_id)
        if reports:
            print export_PSFM(reports[0].get_meta_sites(), format='JASPAR')

def run():
    """Entry point for the script."""
    generate_motif_file()
