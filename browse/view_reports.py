"""This file contains function definitions that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
curation-site-instance objects, they are grouped by TF and species and rendered
properly."""

from django.shortcuts import get_object_or_404
from django.shortcuts import render_to_response
from django.template import RequestContext

import browse.models as models
import browse.motif_report as motif_report

def view_reports_by_TF_and_species(request, TF_id, species_id):
    """Finds sites and generates motif reports given an organism and TF ID."""
    org = get_object_or_404(models.Taxonomy, pk=species_id)
    TF = get_object_or_404(models.TF, TF_id=TF_id)
    
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy=org,
        curation__TF_instances__TF=TF)
    if not request.GET.get('integrate_non_motif', None):
        curation_site_instances = curation_site_instances.exclude(
            site_type='non_motif_associated')

    reports = motif_report.make_distinct_reports(curation_site_instances)
    ensemble_report = motif_report.make_ensemble_report(curation_site_instances)
    return render_to_response(
        'view_reports.html',
        {'reports': [r.generate_view_reports_dict() for r in reports],
         'ensemble_report': ensemble_report.generate_view_reports_dict()},
        context_instance=RequestContext(request))
    
def view_reports_by_id_list(request):
    """Returns motif reports from given Curation_SiteInstance object IDs."""
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        pk__in=request.POST['csi_list'].strip().split(','))

    if request.POST.get('integrate_non_motif', None):
        curation_site_instances = curation_site_instances.exclude(
            site_type='non_motif_associated')

    reports = motif_report.make_distinct_reports(curation_site_instances)
    ensemble_report = motif_report.make_ensemble_report(cur_site_insts)

    return render_to_response(
        'view_reports.html',
        {'reports': [report.generate_view_reports_dict() for report in reports],
         'ensemble_report': ensemble_report.generate_view_reports_dict(),
         'curation_site_instances': ','.join(
             map(lambda csi: '%d'%csi.pk, curation_site_instances)),
         'integrate_non_motif': integrate_non_motif,
         'by_id': True},
        context_instance=RequestContext(request))

def view_reports_by_TF_family(request, object_id): 
    """Returns the motif reports for the given TF family."""
    pass

def view_reports_by_TF(request, object_id):
    """Returns the motif reports for the given TF."""
    pass

def view_reports_by_all_techniques(request, function):
    pass
    
def view_reports_by_technique_category(request, category_function, object_id):
    """Returns the motif reports for a given experimental technique category."""
    pass

def view_reports_by_technique(request, object_id):
    """Returns the motif reports for a given experimental technique."""
    pass

def view_reports_by_taxonomy(request, object_id):
    """Returns the motif reports for a given taxon."""
    pass
