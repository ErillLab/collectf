from collectfapp import models
from collectfapp import bioutils

def get_all_curation_site_instances():
    all_curation_site_instances = list(models.Curation_SiteInstance.objects.all())
    all_curation_site_instances.sort(key=lambda csi: -len(csi.site_instance.seq))
    print ','.join([
        'curation_site_instance_id',
        'curation_id',
        'site_instance_id',
        'site_instance_len',
        'is_motif_associated',
        ])
    for csi in all_curation_site_instances:
        print ','.join(map(str, [
            csi.pk,
            csi.curation.curation_id,
            csi.site_instance,
            len(csi.site_instance.seq),
            csi.is_motif_associated,
            ]))

        
def run():
    get_all_curation_site_instances()
    print "Running misc.py"
    
    
