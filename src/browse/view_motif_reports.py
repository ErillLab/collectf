"""This module contains function that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
Curation_SiteInstance objects, they are grouped by TF and species and rendered
properly.
"""

from django.shortcuts import get_object_or_404
from django.shortcuts import render

from core import models
from .motif_report import build_motif_reports


def render_motif_report(request, curation_site_instances):
    """Renders the motif report page."""
    motif_reports = build_motif_reports(curation_site_instances)
    for motif_report in motif_reports:
        motif_report.meta_sites[0].curation_site_instances
    return render(request, 'view_motif_reports.html',
                  {'motif_reports': motif_reports})


def view_reports_by_TF_and_species(request, TF_id, species_id):
    """Finds sites and generates motif reports given an organism and TF ID."""
    org = get_object_or_404(models.Taxonomy, pk=species_id)
    TF = get_object_or_404(models.TF, TF_id=TF_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy=org,
        curation__TF_instances__TF=TF,
        site_type='motif_associated')
    return render_motif_report(request, curation_site_instances)


def view_reports_by_TF_family(request, object_id):
    """Returns the motif reports for the given TF family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__family=TF_family)
    return render_motif_report(request, curation_site_instances)


def view_reports_by_TF(request, object_id):
    """Returns the motif reports for the given TF."""
    TF = get_object_or_404(models.TF, TF_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF=TF)
    return render_motif_report(request, curation_site_instances)


def view_reports_by_technique_function(request, function):
    """Returns the motif reports for the given technique function."""
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        experimental_techniques__preset_function=function)
    return render_motif_report(request, curation_site_instances)


def view_reports_by_technique_category(request, category_function, object_id):
    """Returns motif reports for a given experimental technique category."""
    category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                 category_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        experimental_techniques__categories=category)
    return render_motif_report(request, curation_site_instances)


def view_reports_by_technique(request, object_id):
    """Returns the motif reports for a given experimental technique."""
    technique = get_object_or_404(models.ExperimentalTechnique,
                                  technique_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        experimental_techniques=technique)
    return render_motif_report(request, curation_site_instances)


def view_reports_by_taxonomy(request, object_id):
    """Returns the motif reports for a given taxon."""
    taxonomy = get_object_or_404(models.Taxonomy, pk=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in=taxonomy.get_all_species())
    return render_motif_report(request, curation_site_instances)


def view_reports_by_curation_site_instance_ids(request):
    """Returns motif reports from given Curation_SiteInstance object IDs."""
    curation_site_instance_ids = [
        x.strip() for x in request.POST['curation_site_instance_list'].split()]
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        pk__in=curation_site_instance_ids)
    return render_motif_report(request, curation_site_instances)
