# This script checks for inconsistency in CollecTF database.
# TODO(sefa): Set up a cron job for this script.

from base import models

def each_curation_has_only_one_tf():
    """Check if each curation have only one TF associated with it."""
    for curation in models.Curation.objects.all():
        tfs = set(tf_instance.TF for
                  tf_instance in curation.TF_instances.all())
        if not len(tfs) == 1:
            print ("Curation %d has %d associated TFs." %
                   (curation.curation_id, len(tfs)))

def tf_instances_have_valid_accession_number():
    """Finds invalid TF instances and curations associated with it, if any."""
    tf_instances = models.TFInstance.objects.all()
    allowed_prefixes = ['NP', 'YP', 'WP']
    for tf_instance in tf_instances:
        if not any(tf_instance.protein_accession.startswith(prefix)
                   for prefix in allowed_prefixes):
            print ("TF instance %s has an invalid accession number." %
                   tf_instance.protein_accession)
            
def run():
    each_curation_has_only_one_tf()
    tf_instances_have_valid_accession_number()
            
