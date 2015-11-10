"""Module for searching binding motifs in CollecTF."""

from django.contrib import messages
from django.core.exceptions import ValidationError
from django.core.urlresolvers import reverse
from django.db.models import Q
from django.http import HttpResponseRedirect
from django.shortcuts import render

from core import models
from .motif_report import build_motif_reports, build_ensemble_report


def search(request):
    """Handler for search view.

    On GET, renders the page with the list of TFs, species and
    experimental techniques.

    On POST, identifies selected TFs, species and techniques and retrives
    Curation_SiteInstance that fits to search criteria.
    """
    if request.POST:
        return search_post(request)

    return search_get(request)


def search_get(request):
    """Renders the page to select TFs, species and techniques."""

    binding, expression, insilico = get_all_techniques()
    return render(request, 'search.html',
                  {'TF_families': get_all_TF_families(),
                   'phyla': get_all_phyla(),
                   'all_techniques': get_all_techniques()})


def get_all_TF_families():
    """Returns all TF families in the database."""
    return models.TFFamily.objects.all().order_by('name')


def get_all_phyla():
    """Returns all phyla in the database."""
    return models.Taxonomy.objects.filter(rank='phylum').order_by('name')


def get_all_techniques():
    """Gets all techniques in the database and group them by type.

    Type is one of {binding, expression, in-silico}.
    """
    all_techniques = {'binding': {}, 'expression': {}, 'insilico': {}}
    for category in models.ExperimentalTechniqueCategory.objects.all():
        techniques = models.ExperimentalTechnique.objects.filter(
            categories=category).order_by('name')
        binding_techniques = techniques.filter(preset_function='binding')
        if binding_techniques:
            all_techniques['binding'][category.name] = binding_techniques
        expression_techniques = techniques.filter(preset_function='expression')
        if expression_techniques:
            all_techniques['expression'][category.name] = expression_techniques
        insilico_techniques = techniques.filter(preset_function='insilico')
        if insilico_techniques:
            all_techniques['insilico'][category.name] = insilico_techniques
    return all_techniques


def search_post_helper(request):
    """Handler for search_post requests.

    Given request, gets motif- and non-motif-associated sites from database and
    renders the results.
    """

    def get_TF_input():
        """Gets the TF input from the form."""
        TF_input = [TF_id for TF_id in request.POST.getlist('tf_input')
                    if TF_id != 'on']
        if not TF_input:
            raise ValidationError("Select at least one TF or TF family.")
        return models.TF.objects.filter(TF_id__in=TF_input)

    def get_species_input():
        """Gets the species input from the form."""
        species_input = [pk for pk in request.POST.getlist('species_input')
                         if pk != 'on']
        if not species_input:
            raise ValidationError("Select at least one taxonomic group.")
        return models.Taxonomy.objects.filter(pk__in=species_input)

    def get_technique_input():
        """Gets the experimental technique input from the form."""
        # Get category inputs
        category_inputs = [
            [cat_id for cat_id in request.POST.getlist(field) if cat_id != 'on']
            for field in ['cat_input_1', 'cat_input_2', 'cat_input_3']]
        techniques = [models.ExperimentalTechnique.objects.filter(
            pk__in=category_input) for category_input in category_inputs]
        if not any(techniques):
            raise ValidationError("Select at least one experimental technique.")
        return techniques

    def techniques_to_Q(techniques):
        """Creates a query filter for the given techniques.

        Returns a query filter for objects that are associated with any of the
        given techniques.
        """
        q = Q(curation__curation_id=-9999)
        for t in techniques:
            q = q | Q(experimental_techniques=t)
        return q

    def filter_by_technique(curation_site_instances):
        """Filters Curation_SiteInstance objects."""
        boolean1 = request.POST['boolean1']
        boolean2 = request.POST['boolean2']
        # Create filter queries
        q1, q2, q3 = map(techniques_to_Q, get_technique_input())
        if boolean1 == 'and' and boolean2 == 'and':
            curation_site_instances = curation_site_instances.filter(q1, q2, q3)
        elif boolean1 == 'and' and boolean2 == 'or':
            # (A and B) or C <-> (A or C) and (B or C)
            curation_site_instances = curation_site_instances.filter(q1|q2, q2|q3)
        elif boolean1 == 'or' and boolean2 == 'and':
            curation_site_instances = curation_site_instances.filter(q1|q2).filter(q3)
        elif boolean1 == 'or' and boolean2 == 'or':
            curation_site_instances = curation_site_instances.filter(q1|q2|q3)
        else:
            assert False, 'Shouldnt be here, unhandled case.'
        return curation_site_instances

    # Get inputs
    TFs = get_TF_input()
    species = get_species_input()
    # Get curation-site-instance objects, filtered by TF and species
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__in=TFs,
        site_instance__genome__taxonomy__in=species)
    # Filter Curation-site-instance objects by experimental techniques
    curation_site_instances = filter_by_technique(curation_site_instances)
    return curation_site_instances.distinct()


def search_post(request):
    """Handler for database search request."""
    try:
        curation_site_instances = search_post_helper(request)
        motif_reports = build_motif_reports(curation_site_instances)
        motif_ensemble_report = build_ensemble_report(curation_site_instances)
        return render(
            request, 'search_results.html',
            {'motif_reports': motif_reports,
             'ensemble_motif_report': motif_ensemble_report})

    except ValidationError, e:
        messages.add_message(request, messages.ERROR, e.message)
        return HttpResponseRedirect(reverse(search))
