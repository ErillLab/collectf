from browse_base import *

def browse_by_site(request, dbxref_id):
    """Handler for browsing site instances."""
    site_instance_id = utils.dbxref2id(dbxref_id)
    site_instance = models.SiteInstance.objects.get(site_id=site_instance_id)
    # get curations related to this site instance
    curations = models.Curation.objects.filter(site_instances=site_instance)
    # get curation_site_instances related to this site instance (needed for regulation)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(site_instance=site_instance)
    
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
                  {'site_instance': site_instance,
                   'curations': curations,
                   'regulations': regulations,
                   'dbxref': utils.id2dbxref(site_instance.site_id)
                  },
                  context_instance=RequestContext(request))
