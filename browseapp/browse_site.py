from django.contrib import messages
from browse_base import *
import django.views.defaults


def browse_by_site(request, dbxref_id):
    """Handler for browsing site instances."""
    id = utils.dbxref2id(dbxref_id)
    curation_site_instance = models.Curation_SiteInstance.objects.get(pk=id)

    # get all curation-site-instance objects overlaps with this one
    all_curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instance=curation_site_instance.curation.TF_instance,
        site_instance__genome=curation_site_instance.site_instance.genome,
        is_motif_associated=True)

    curation_site_instances = [csi for csi in all_curation_site_instances
                               if bioutils.overlap_site_meta_site(curation_site_instance, [csi])]
    
    # get curations
    curation_ids = list(set(csi.curation.curation_id for csi in curation_site_instances))
    curations = models.Curation.objects.filter(curation_id__in=curation_ids)
    # get regulations
    regulations = []
    for curation_site_instance in curation_site_instances:
        for reg in curation_site_instance.regulation_set.all():
            for r in regulations: # check if the same gene is already in the list
                if r.gene == reg.gene:
                    if r.evidence_type=='inferred' and reg.evidence_type=='exp_verified':
                        r.evidence_type='exp_verified'
                    break
            else:
                regulations.append(reg)
                
    return render(request,
                  "browse_site.html",
                  {
                      'head_curation_site_instance': curation_site_instance,
                      'curation_site_instances': curation_site_instances,
                      'curations': curations,
                      'regulations': regulations,
                      'dbxref': utils.id2dbxref(curation_site_instance.pk)
                  },
                  context_instance=RequestContext(request))
