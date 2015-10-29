"""This module contains function that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
Curation_SiteInstance objects, they are grouped by TF and species and rendered
properly.
"""

from django.shortcuts import get_object_or_404
from django.shortcuts import redirect
from django.shortcuts import render_to_response
from django.template import RequestContext

from browse import models
from browse.motif_report import make_reports
from browse.motif_report import make_ensemble_report
from browse.dbxref import from_uniprot_dbxref
    
def view_reports_by_id_list(request):
    """Returns motif reports from given Curation_SiteInstance object IDs."""
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        pk__in=request.POST['csi_list'].strip().split(','))
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)

    return render_to_response(
        'view_reports.html',
        {'reports': reports,
         'ensemble_report': ensemble_report,
         'by_id': True},
        context_instance=RequestContext(request))

def render_report_to_response(request, reports, ensemble_report):
    """Helps rendering MotifReport object to the response."""
    return render_to_response(
        'view_reports.html',
        {'reports': reports, 'ensemble_report': ensemble_report},
        context_instance=RequestContext(request))

def view_reports_by_uniprot_dbxref(request, uniprot_dbxref):
    """Builds motif report for the given UniProt dbxref identifier."""
    TF_instance_id = from_uniprot_dbxref(uniprot_dbxref)
    return redirect(view_reports_by_TF_instance, TF_instance_id)

def view_reports_by_uniprot_accession(request, uniprot_accession):
    """Builds motif report for the given UniProt accession number."""
    TF_instance = get_object_or_404(models.TFInstance,
                                    uniprot_accession=uniprot_accession)
    return redirect(view_reports_by_TF_instance, TF_instance.TF_instance_id)

def view_reports_by_TF_instance(request, TF_instance_id):
    """Finds sites and generates motif reports given a TF instance."""
    TF_instance = get_object_or_404(models.TFInstance,
                                    TF_instance_id=TF_instance_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances=TF_instance,
        site_type='motif_associated')
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)
    
def view_reports_by_TF_and_species(request, TF_id, species_id):
    """Finds sites and generates motif reports given an organism and TF ID."""
    org = get_object_or_404(models.Taxonomy, pk=species_id)
    TF = get_object_or_404(models.TF, TF_id=TF_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy=org,
        curation__TF_instances__TF=TF,
        site_type='motif_associated')
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_TF_family(request, object_id): 
    """Returns the motif reports for the given TF family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__family=TF_family)
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_TF(request, object_id):
    """Returns the motif reports for the given TF."""
    TF = get_object_or_404(models.TF, TF_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF=TF)
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_all_techniques(request, function):
    """Returns the motif reports for the given technique function."""
    techniques = models.ExperimentalTechnique.objects.filter(
        preset_function=function)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        experimental_techniques=techniques)
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)
    
def view_reports_by_technique_category(request, category_function, object_id):
    """Returns the motif reports for a given experimental technique category."""
    category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                 category_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        experimental_techniques_category=category)
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_technique(request, object_id):
    """Returns the motif reports for a given experimental technique."""
    technique = get_object_or_404(models.ExperimentalTechnique,
                                  technique_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        experimental_techniques=technique)
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)

def view_reports_by_taxonomy(request, object_id):
    """Returns the motif reports for a given taxon."""
    taxonomy = get_object_or_404(models.Taxonomy, pk=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in=taxonomy.get_all_species())
    reports = make_reports(curation_site_instances)
    ensemble_report = make_ensemble_report(curation_site_instances)
    return render_report_to_response(request, reports, ensemble_report)
