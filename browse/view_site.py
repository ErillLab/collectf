from django.contrib import messages
from base import models
from base import metasite
from base import bioutils
from base.templatetags import dbxref
from django.shortcuts import render
from django.template import RequestContext

def view_site(request, dbxref_id):
    """Handler to view site instances."""
    try:
        id = dbxref.dbxref2id(dbxref_id)
        curation_site_instance = models.Curation_SiteInstance.objects.get(pk=id)
        if curation_site_instance.is_obsolete:
            messages.add_message(request, messages.ERROR,
                                 "The requested site instance seems obsolete. "
                                 "It will be removed on the next release. "
                                 "Description: %s" % curation_site_instance.why_obsolete)
            
    except:
        messages.add_message(request, messages.ERROR, "Invalid binding site id.")
        return HttpResponseRedirect(reverse(base.views.home))

    # get all curation-site-instance objects overlaps with this one
    all_curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances=curation_site_instance.curation.TF_instances,
        site_instance__genome=curation_site_instance.site_instance.genome,
        site_type=["motif_associated", "var_motif_associated"]).all()

    # get all curation_site_instances that overlap with the queried site
    meta_site = metasite.MetaSite(curation_site_instance)
    for csi in all_curation_site_instances:
        if csi==curation_site_instance: continue
        if meta_site.membership_test(csi):
                meta_site.add_cur_site_inst(csi)

    alignment = None
    if len(meta_site.cur_site_insts) > 1:
        alignment = bioutils.run_lasagna(map(lambda csi: csi.site_instance, meta_site.cur_site_insts), trim=False)
    
    return render(request,
                  "view_site.html",
                  {
                      'head_csi': curation_site_instance,
                      'csis': meta_site.cur_site_insts,
                      'regulations': meta_site.regulations,
                      'dbxref': dbxref.id2dbxref(curation_site_instance.pk),
                      'alignment': alignment,
                  },
                  context_instance=RequestContext(request))
