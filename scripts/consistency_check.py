# This script checks for inconsistency in CollecTF database.
# TODO(sefa): Set up a cron job for this script.

from base import models

def each_curation_has_only_one_tf():
    """Check if each curation have only one TF associated with it."""
    for curation in models.Curation.objects.all():
        tfs = set(tf_instance.TF for
                  tf_instance in curation.TF_instances.all())
        if not len(tfs) == 1:
            print ("Curation %d has more than %d associated TFs." %
                   (curation.curation_id, len(tfs)))

def run():
    each_curation_has_only_one_tf()
            
