"""This module contains function that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
Curation_SiteInstance objects, they are grouped by TF and species and rendered
properly.
"""

from django.shortcuts import get_object_or_404
from django.shortcuts import redirect
from django.shortcuts import render_to_response
from django.template import RequestContext

import browse.models as models
import browse.motif_report as motif_report
from browse.static_reports import get_static_reports
from browse.dbxref import from_uniprot_dbxref
    
def view_reports_by_id_list(request):
    """Returns motif reports from given Curation_SiteInstance object IDs."""
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        pk__in=request.POST['csi_list'].strip().split(','))
    reports = motif_report.make_distinct_reports(curation_site_instances)
    ensemble_report = motif_report.make_ensemble_report(cur_site_insts)

    return render_to_response(
        'view_reports.html',
        {'reports': [report.generate_view_reports_dict() for report in reports],
         'ensemble_report': ensemble_report.generate_view_reports_dict(),
         'curation_site_instances': ','.join(
             map(lambda csi: '%d'%csi.pk, curation_site_instances)),
         'by_id': True},
        context_instance=RequestContext(request))

def render_report_to_response(request, reports, ensemble_report):
    """Helps rendering MotifReport object to the response."""
    return render_to_response(
        'view_reports.html',
        {'reports': [r.generate_view_reports_dict() for r in reports],
         'ensemble_report': ensemble_report.generate_view_reports_dict()},
        context_instance=RequestContext(request))

def view_reports_by_uniprot_id(request, uniprot_dbxref):
    """Builds motif report for the given UniProt dbxref identifier."""
    TF_instance_id = from_uniprot_dbxref(uniprot_dbxref)
    return redirect(view_reports_by_TF_instance, TF_instance_id)

def view_reports_by_TF_instance(request, TF_instance_id):
    """Finds sites and generates motif reports given a TF instance."""
    TF_instance = get_object_or_404(models.TFInstance,
                                    TF_instance_id=TF_instance_id)
    reports, ensemble_report = get_static_reports(
        'tf_instance_%s' % TF_instance_id)
    return render_report_to_response(request, reports, ensemble_report)
    
def view_reports_by_TF_and_species(request, TF_id, species_id):
    """Finds sites and generates motif reports given an organism and TF ID."""
    org = get_object_or_404(models.Taxonomy, pk=species_id)
    TF = get_object_or_404(models.TF, TF_id=TF_id)
    reports, ensemble_report = get_static_reports(
        'tf_%s_species_%s' % (TF_id, species_id))
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_TF_family(request, object_id): 
    """Returns the motif reports for the given TF family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    reports, ensemble_report = get_static_reports('TF_family_%s' % object_id)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_TF(request, object_id):
    """Returns the motif reports for the given TF."""
    TF = get_object_or_404(models.TF, TF_id=object_id)
    reports, ensemble_report = get_static_reports('TF_%s' % object_id)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_all_techniques(request, function):
    """Returns the motif reports for the given technique function."""
    reports, ensemble_report = get_static_reports(
        'experimental_technique_all_' + function)
    return render_report_to_response(request, reports, ensemble_report)
    
def view_reports_by_technique_category(request, category_function, object_id):
    """Returns the motif reports for a given experimental technique category."""
    category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                 category_id=object_id)
    reports, ensemble_report = get_static_reports(
        'experimental_technique_category_%s' % object_id)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_technique(request, object_id):
    """Returns the motif reports for a given experimental technique."""
    technique = get_object_or_404(models.ExperimentalTechnique,
                                  technique_id=object_id)
    reports, ensemble_report = get_static_reports(
        'experimental_technique_%s' % object_id)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_taxonomy(request, object_id):
    """Returns the motif reports for a given taxon."""
    taxonomy = get_object_or_404(models.Taxonomy, taxononmy_id=object_id)
    reports, ensemble_report = get_static_reports('taxonomy_%s' % object_id)
    return render_report_to_response(request, reports, ensemble_report)
