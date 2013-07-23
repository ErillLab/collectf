from django.contrib import messages
from browse_base import *
import django.views.defaults


def browse_by_site(request, dbxref_id):
    """Handler for browsing site instances."""
    id = utils.dbxref2id(dbxref_id)
    meta_site_instance = models.MetaSiteInstance.objects.get(pk=id)
    
    # get curation_site_instance related to this metasiteinstance
    curation_site_instances = models.Curation_SiteInstance.objects.filter(meta_site_instance=meta_site_instance).distinct()
    # get curations
    curation_ids = curation_site_instances.values_list('curation', flat=True)
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
                      'meta_site_instance': meta_site_instance,
                      'curation_site_instances': curation_site_instances,
                      'curations': curations,
                      'regulations': regulations,
                      'dbxref': utils.id2dbxref(meta_site_instance.pk)
                  },
                  context_instance=RequestContext(request))
