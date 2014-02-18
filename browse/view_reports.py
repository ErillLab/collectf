"""This file contains function definitions that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
curation-site-instance objects, they are grouped by TF and species and rendered
properly."""

from django.shortcuts import render
from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.template import RequestContext
from django.core.urlresolvers import reverse
import base
import models
import motif_report

def view_reports(request):
    """Handler for the post request, performed on browse by TF/species/technique
    page. Given the list of motif-associated and non-motif-associated
    curation-site-instance ids, they are grouped by TF and species, and
    rendered"""
    # make sure that it is a POST request
    if not request.POST: return HttpResponseRedirect(reverse(base.views.home))

    cur_site_inst_ids = request.POST["csi_list"].strip().split(',')

    # integrate non-motif: if True, the experimental evidence and regulation
    # information from non-motif associated sites are integrated into the report
    integrate_non_motif = bool("integrate_non_motif" in request.POST)

    # get motif-associated and non-motif-associated curation-site-instances and make reports
    cur_site_insts = models.Curation_SiteInstance.objects.filter(pk__in=cur_site_inst_ids)
    reports = motif_report.make_reports(cur_site_insts)
    ensemble_report = motif_report.make_ensemble_report(cur_site_insts)
    return render_to_response("view_reports.html",
                              {
                                  "reports": [report.generate_view_reports_dict()
                                              for report in reports],
                                  "ensemble_report": ensemble_report.generate_view_reports_dict(),
                                  "cur_site_insts": ','.join(map(lambda x: str(x.pk), cur_site_insts))
                              },
                              context_instance=RequestContext(request))


    
    
    
