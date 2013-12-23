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

def get_all_curations_with_chip():
    all_curations = models.Curation.objects.filter(experimental_techniques__name__in=['ChIP-chip','ChIP-Seq','ChIP-PCR']).distinct()
    print 'curation_id, PMID, used_techniques'
    for curation in all_curations:
        print curation.curation_id, '\t',
        print curation.publication.pmid, '\t',
        print map(str, curation.experimental_techniques.values_list('name', flat=True))

def get_all_TF_instances_without_curation():
    used_TF_instances = models.Curation.objects.values_list('TF_instance',flat=True).distinct()
    not_used = models.TFInstance.objects.exclude(protein_accession__in=used_TF_instances)
    for x in not_used:
        print x.protein_accession,

def foo():
    num_all = models.Curation_SiteInstance.objects.all().count()
    ma = models.Curation_SiteInstance.objects.filter(is_motif_associated=True)
    nma = models.Curation_SiteInstance.objects.filter(is_motif_associated=False)
    print num_all, ma.count(), nma.count()
    for csi in ma:
        csi.site_type = "motif_associated"
        csi.save()
    for csi in nma:
        csi.site_type = "non_motif_associated"
        csi.save()

    

        
def run():
    #get_all_curation_site_instances()
    #get_all_curations_with_chip()
    #get_all_TF_instances_without_curation()
    foo()
    
    
