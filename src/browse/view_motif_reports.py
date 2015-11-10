"""This module contains function that are used to display generated
reports. Given a collection of motif-associated and non-motif-associated
Curation_SiteInstance objects, they are grouped by TF and species and rendered
properly.
"""

from django.shortcuts import get_object_or_404
from django.shortcuts import get_list_or_404
from django.shortcuts import render

from core import models
from . import motif_report


def group_curation_site_instances(curation_site_instances):
    """Groups Curation_SiteInstance objects by TF-instance and genome."""
    TF_genome_pairs = curation_site_instances.filter(
        site_type='motif_associated').values_list(
        'site_instance__genome', 'curation__TF_instances').distinct()
    return [curation_site_instances.filter(site_instance__genome=genome,
                                           curation__TF_instances=TF_instances)
            for genome, TF_instances in TF_genome_pairs]


def render_motif_report(request, curation_site_instances):
    """Renders the motif report page."""
    motif_reports = build_motif_reports(curation_site_instances)
    for motif_report in motif_reports:
        motif_report.meta_sites[0].techniques
    return render(request, 'view_motif_reports.html',
                  {'motif_reports': motif_reports})


def build_motif_reports(curation_site_instances):
    """Given a list of Curation_SiteInstance objects, returns MotifReports."""
    curation_site_instance_grps = group_curation_site_instances(
        curation_site_instances)
    for grp in curation_site_instance_grps:
        print len(grp)
    return [motif_report.MotifReport(curation_site_instance_grp)
            for curation_site_instance_grp in curation_site_instance_grps]


def build_ensemble_report(curation_site_instances):
    """Given Curation_SiteInstance objects, returns the EnsembleMotifReport."""
    return motif_report.EnsembleMotifReport(curation_site_instances)


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
    """Returns the motif reports for a given experimental technique category."""
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
