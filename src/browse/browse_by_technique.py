"""The view functions for browsing by TF and TF family."""

from django.core.urlresolvers import reverse
from django.shortcuts import render
from django.shortcuts import get_object_or_404

from core import models
from .browse_by_utils import curation_site_instances_values_list


def browse_technique(request):
    """Returns all techniques, grouped by categories."""
    # Group all binding and expression techniques by their category.
    categories = models.ExperimentalTechniqueCategory.objects.all()
    techniques_dict = {'binding': {}, 'expression': {}}
    for category in categories:
        techniques = models.ExperimentalTechnique.objects.filter(
            categories=category)
        binding_techs = techniques.filter(preset_function='binding')
        if binding_techs:
            techniques_dict['binding'][category] = binding_techs
        expression_techs = techniques.filter(preset_function='expression')
        if expression_techs:
            techniques_dict['expression'][category] = expression_techs
    return render(request, 'browse_by_technique.html',
                  {'techniques_dict': techniques_dict})


def get_results_by_technique_function(request, function):
    """Returns motif reports of sites validated by binding or expression."""
    assert function in ['binding', 'expression']

    techniques = models.ExperimentalTechnique.objects.filter(
        preset_function=function)
    TF_species_pairs = curation_site_instances_values_list(
        models.Curation_SiteInstance.objects.filter(
            experimental_techniques__in=techniques))
    title_lookup = {'binding': 'Detection of binding',
                    'expression': 'Assessment of expression'}
    return render(request, 'browse_results.html',
                  {'title': title_lookup[function],
                   'description': '',
                   'TF_species_pairs': TF_species_pairs,
                   'ensemble_report_url': reverse(
                       'view_motif_reports_by_technique_function',
                       args=(function,))})


def get_results_by_technique_category(request, category_function, object_id):
    """Returns motif reports by binding/expression category."""
    assert category_function in ['binding', 'expression']
    category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                 category_id=object_id)
    techniques = models.ExperimentalTechnique.objects.filter(
        categories=category, preset_function=category_function)
    TF_species_pairs = curation_site_instances_values_list(
        models.Curation_SiteInstance.objects.filter(
            experimental_techniques__in=techniques))
    return render(request, 'browse_results.html',
                  {'title': category.name,
                   'description': category.description,
                   'TF_species_pairs': TF_species_pairs,
                   'ensemble_report_url': reverse(
                       'view_motif_reports_by_technique_category',
                       args=(category_function, object_id))})


def get_results_by_technique(request, object_id):
    """Returns motif reports by experimental technique ID."""
    technique = get_object_or_404(models.ExperimentalTechnique,
                                  technique_id=object_id)
    TF_species_pairs = curation_site_instances_values_list(
        models.Curation_SiteInstance.objects.filter(
            experimental_techniques=technique))
    return render(request, 'browse_results.html',
                  {'title': technique.name,
                   'description': technique.description,
                   'TF_species_pairs': TF_species_pairs,
                   'ensemble_report_url': reverse(
                       'view_motif_reports_by_technique', args=(object_id,))})
