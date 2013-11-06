from django.contrib import messages
from browse_base import *
import django.views.defaults

import collectfapp.views

def browse_by_site(request, dbxref_id):
    """Handler for browsing site instances."""
    try:
        id = dbxref_utils.dbxref2id(dbxref_id)
        curation_site_instance = models.Curation_SiteInstance.objects.get(pk=id)
        if curation_site_instance.is_obsolete:
            messages.add_message(request, messages.ERROR,
                                 "The requested site instance seems obsolete. "
                                 "It will be removed on the next release. "
                                 "Description: %s" % curation_site_instance.why_obsolete)
            
    except:
        messages.add_message(request, messages.ERROR, "Invalid binding site id.")
        return HttpResponseRedirect(reverse(collectfapp.views.home))

    # get all curation-site-instance objects overlaps with this one
    all_curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instance=curation_site_instance.curation.TF_instance,
        site_instance__genome=curation_site_instance.site_instance.genome,
        is_motif_associated=True).all()

    # get all curation_site_instances that overlap with the queried site
    curation_site_instances = [curation_site_instance]
    for csi in all_curation_site_instances:
        if csi==curation_site_instance: continue
        if bioutils.overlap_site_meta_site(csi, curation_site_instances):
            curation_site_instances.append(csi)


    alignment = None
    if len(curation_site_instances) > 1:
        alignment = bioutils.call_lasagna(map(lambda csi: csi.site_instance, curation_site_instances), trim=False)
    
    # get curations
    curation_ids = models.Curation.objects.values_list('curation_id', flat=True).distinct()
    curations = models.Curation.objects.filter(curation_id__in=curation_ids)
    # get regulations
    regulations = []
    for csi in curation_site_instances:
        for reg in csi.regulation_set.all():
            for r in regulations: # check if the same gene is already in the list
                if r.gene == reg.gene:
                    if r.evidence_type=='inferred' and reg.evidence_type=='exp_verified':
                        r.evidence_type='exp_verified'
                    break
            else:
                regulations.append(reg)
                
    return render(request,
                  "view_site.html",
                  {
                      'head_csi': curation_site_instance,
                      'csis': curation_site_instances,
                      'regulations': regulations,
                      'dbxref': dbxref_utils.id2dbxref(curation_site_instance.pk),
                      'alignment': alignment,
                  },
                  context_instance=RequestContext(request))
