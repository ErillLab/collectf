"""
This script checks for inconsistency in CollecTF database.
"""
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
    """Checks if all TF-instances have valid accession numbers."""
    tf_instances = models.TFInstance.objects.all()
    allowed_prefixes = ['NP', 'YP', 'WP']
    for tf_instance in tf_instances:
        if not any(tf_instance.protein_accession.startswith(prefix)
                   for prefix in allowed_prefixes):
            print ("TF instance %s has an invalid accession number." %
                   tf_instance.protein_accession)

def curations_with_no_sites():
    """Checks if there are curations with no sites."""
    for curation in models.Curation.objects.all():
        if not curation.site_instances.all():
            print "Curation %d has no site instances." % curation.curation_id

def curation_site_instances_link_to_curation():
    """Checks if all Curation_SiteInstance objects link to Curation objects."""
    for curation_site_instance in models.Curation_SiteInstance.objects.all():
        if not curation_site_instance.curation:
            print ("Curation_SiteInstance %d doesn't have link to a curation." %
                   curation_site_instance.pk)
            print ("Associated curations are:",
                   models.Curation.objects.filter(TF_instances=tf_instance))

def tf_instances_have_associated_tf():
    """Checks if all TF-instances have an associated TF."""
    for tf_instance in models.TFInstance.objects.all():
        if not tf_instance.TF:
            print ("TF instance %s doesn't have a TF." %
                   tf_instance.uniprot_accession)

def run():
    """Entry point for the script."""
    each_curation_has_only_one_tf()
    tf_instances_have_valid_accession_number()
    curations_with_no_sites()
    curation_site_instances_link_to_curation()
    tf_instances_have_associated_tf()
