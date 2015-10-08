"""
Scripts for CollecTF database schema changes. Intented to run once.
"""

from tqdm import tqdm

from base import models

def move_TF_from_curation_to_TF_instance_table():
    """Moves relationship between curation and TF to TF and TF-instance."""
    for curation in tqdm(models.Curation.objects.all()):
        for TF_instance in curation.TF_instances.all():
            TF_instance.TF = curation.TF
            TF_instance.save()
    
def run():
    """Entry point for the script."""
    move_TF_from_curation_to_TF_instance_table()
