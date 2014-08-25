# This script file contains some functions that were used to make some changes
# (moving the data from one table to the other, etc.) in the database scheme
# without losing any data.

# They are mostly run only once, but are kept here for future reference.

from base import models

def update_db():
    """Apply following changes
    - move TF-function from curation to curation-site-instance
    - change curation->TF-instance field to many2many field
    - move experimental-techniques from curation to curation-site-instance
    """
    cur_site_insts = models.Curation_SiteInstance.objects.all()
    for cur_site_inst in cur_site_insts:
        curation = cur_site_inst.curation
        print curation
        # Update TF function
        cur_site_inst.TF_function = curation.TF_function
        # Update experimental techniques
        for exp in curation.experimental_techniques.all():
           cur_site_inst.experimental_techniques.add(exp)
        # Update TF type
        cur_site_inst.TF_type = curation.TF_type
        cur_site_inst.save()
    print "all curation-site-intances updated."

    
    # Update curation TF instances
    all_curations = models.Curation.objects.all()
    for cur in all_curations:
        cur.TF_instances.add(cur.TF_instance)
        cur.save()
    print "all curations updated."

def update_motif_ids():
    for cur_site_inst in models.Curation_SiteInstance.objects.all():
        if cur_site_inst.site_type == 'motif_associated':
            cur_site_inst.motif_id = 0
        else:
            cur_site_inst.motif_id = -1

        cur_site_inst.save()


def run():
    update_motif_ids()
