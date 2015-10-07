"""
The view functions for browsing by TF and TF family.
"""

from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404, get_list_or_404
from django.shortcuts import render_to_response
from django.template import RequestContext

from browse import models
from browse.motif_report import make_reports
from browse.view_reports import view_reports_by_TF_family

def browse_TF(request):
    """Returns the TF tree-view."""
    TF_families = models.TFFamily.objects.all().order_by('name')
    return render_to_response('browse_TF.html',
                              {'TF_families': TF_families},
                              context_instance=RequestContext(request))

def get_tf_instance(accession):
    return get_object_or_404(models.TFInstance, protein_accession=accession)

def get_results_TF_family(request, object_id):
    """Returns TF-species pairs that have binding sites for the given family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    TFs = get_list_or_404(models.TF, family=TF_family)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__in=TFs)
    reports = make_reports(curation_site_instances)
    return render_to_response(
        'browse_results.html',
        {'title': TF_family.name,
         'description': TF_family.description,
         'reports': [r.generate_browse_result_dict() for r in reports],
         'combined_report_url': reverse(view_reports_by_TF_family,
                                        args=(TF_family.TF_family_id,))},
        context_instance=RequestContext(request))

def get_results_TF(request, object_id):
    """Returns TF-species pairs that have binding sites for the given TF."""
    TF = get_object_or_404(models.TF, TF_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__in=TFs)
    reports = make_reports(curation_site_instances)
    return render_to_response(
        'browse_results.html',
        {'title': TF.name,
         'description': TF.description,
         'reports': [r.generate_browse_result_dict() for r in reports]},
        context_instance=RequestContext(request))

