"""This file contains the view functions for browsing by TF"""

from django.shortcuts import render
from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.template import RequestContext
import models
import Queue
import motif_report

def browse_TF(request):
    """View function for browse by TF. Returns the TF tree"""
    TF_families = models.TFFamily.objects.all().order_by('name')
    return render_to_response('browse_TF.html', {'TF_families': TF_families},
                              context_instance=RequestContext(request))

def get_results_TF(request, type_, id):
    """GIven the type (TF or TF family) and the id of the object, return query
    results and list TF/species that have binding site data for the selected TF
    or TF family."""
    if type_ == 'TF':
        TFs = [models.TF.objects.get(TF_id=id)]
        name = TFs[0].name
        desc = TFs[0].description
    elif type_ == 'TF_family':
        # Get all TFs under that family
        TF_family = models.TFFamily.objects.get(TF_family_id=id)
        TFs = models.TF.objects.filter(family=TF_family).all()
        name = TF_family.name
        desc = TF_family.description

    # get all curation-site-instance objects for browsed TFs
    assert TFs
    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        curation__TF__in=TFs)

    # generate all reports
    reports = motif_report.make_reports(cur_site_insts)

    return render_to_response("browse_results.html",
                              {'title': name,
                               'description': desc,
                               'all_cur_site_insts': [pk for report in reports
                                                      for pk in report.get_all_cur_site_insts()],
                               'reports': [report.generate_browse_result_dict() for report in reports],},
                                context_instance=RequestContext(request))
