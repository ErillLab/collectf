"""This file contains views for search function."""

from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages
from django.template import RequestContext
from django.db.models import Q
import browse.models as models
import browse.motif_report as motif_report

def search(request):
    """Handler for search view.

    If it is GET request, renders the page with the list of TFs, species and
    experimental techniques. If it is POST request, identifies selected TFs,
    species and techniques and retrives Curation_SiteInstance that fits to
    search criteria.
    """
    return search_get(request) if not request.POST else search_post(request)

def search_get(request):
    """Renders the page to select TFs, species and techniques."""
    binding, expression, insilico = get_all_techniques()
    return render_to_response('search_.html',
                              {'TF_families': get_TF_families(),
                               'phyla': get_all_phyla(),
                               'binding_techniques': binding,
                               'expression_techniques': expression,
                               'insilico_techniques': insilico},
                              context_instance=RequestContext(request))

def get_TF_families():
    """Returns all TF families in the database."""
    return models.TFFamily.objects.all().order_by('name')

def get_all_phyla():
    """Returns all phyla in the database."""
    return models.Taxonomy.objects.filter(rank='phylum').order_by('name')

def get_all_techniques():
    """Gets all techniques in the database and group them by type.
    
    type is one of {binding, expression, in-silico}
    """
    binding_techniques = {}
    expression_techniques = {}
    insilico_techniques = {}
    all_categories = models.ExperimentalTechniqueCategory.objects.all()
    for category in all_categories:
        techs = models.ExperimentalTechnique.objects.filter(
            categories=category).order_by('name')
        binding_techniques[category.name] = techs.filter(
            preset_function='binding')
        expression_techniques[category.name] = techs.filter(
            preset_function='expression')
        insilico_techniques[category.name] = techs.filter(
            preset_function='insilico')
    # Remove empty keys
    remove_empty_keys = lambda d: dict((x, y) for x, y in d.items() if y)
    binding_techniques = remove_empty_keys(binding_techniques)
    expression_techniques = remove_empty_keys(expression_techniques)
    insilico_techniques = remove_empty_keys(insilico_techniques)
    return binding_techniques, expression_techniques, insilico_techniques

def search_post_helper(request):
    """Handler for search_post requests.

    Given request, get motif- and non-motif-associated sites from database and
    render the results.
    """
    def get_TF_input():
        TF_input = [x for x in request.POST.getlist('tf_input') if x != 'on']
        if not TF_input:
            raise RuntimeError()
        return models.TF.objects.filter(pk__in=TF_input)

    def get_species_input():
        species_input = [x for x in request.POST.getlist('species_input')
                         if x != 'on']
        if not species_input:
            raise RuntimeError()
        return models.Taxonomy.objects.filter(pk__in=species_input)

    def get_technique_input():
        # Get category inputs
        cat_input_1 = [x for x in request.POST.getlist('cat_input_1')
                       if x != 'on']
        cat_input_2 = [x for x in request.POST.getlist('cat_input_2')
                       if x != 'on']
        cat_input_3 = [x for x in request.POST.getlist('cat_input_3')
                       if x != 'on']
        techniques1 = models.ExperimentalTechnique.objects.filter(
            pk__in=cat_input_1)
        techniques2 = models.ExperimentalTechnique.objects.filter(
            pk__in=cat_input_2)
        techniques3 = models.ExperimentalTechnique.objects.filter(
            pk__in=cat_input_3)
        if not (techniques1 or techniques2 or techniques3):
            raise RuntimeError()
        return (techniques1, techniques2, techniques3)

    def techniques_to_Q(techniques):
        """Creates a query filter for the given techniques.

        Returns a query filter for objects that are associated with any of the
        given techniques.
        """
        q = Q(curation__curation_id=-9999)
        for t in techniques:
            q = q | Q(experimental_techniques=t)
        return q

    def filter_by_technique(cur_site_insts):
        """Filters Curation_SiteInstance objects."""
        boolean1 = request.POST['boolean1']
        boolean2 = request.POST['boolean2']
        # Create filter queries
        q1, q2, q3 = map(techniques_to_Q, get_technique_input())
        if boolean1 == 'and' and boolean2 == 'and':
            cur_site_insts = cur_site_insts.filter(q1, q2, q3)
        elif boolean1 == 'and' and boolean2 == 'or':
            # (A and B) or C <-> (A or C) and (B or C)
            cur_site_insts = cur_site_insts.filter(q1|q2, q2|q3)
        elif boolean1 == 'or' and boolean2 == 'and':
            cur_site_insts = cur_site_insts.filter(q1|q2).filter(q3)
        elif boolean1 == 'or' and boolean2 == 'or':
            cur_site_insts = cur_site_insts.filter(q1|q2|q3)
        else:
            assert False, 'Shouldnt be here, unhandled case.'
        return cur_site_insts

    # Get inputs
    TFs = get_TF_input()
    species = get_species_input()
    # Get curation-site-instance objects, filtered by TF and species
    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__in=TFs,
        site_instance__genome__taxonomy__in=species)
    # Filter Curation-site-instance objects by experimental techniques
    cur_site_insts = filter_by_technique(cur_site_insts)
    return cur_site_insts

def search_post(request):
    """Handler for database search request."""
    def raise_validation_error(msg):
        """Checks if all three steps have been submitted."""
        messages.add_message(request, messages.ERROR, message)
        return HttpResponseRedirect(reverse(search))

    try:
        cur_site_insts = search_post_helper(request)
        # make reports
        reports = motif_report.make_reports(cur_site_insts)
        # make ensemble reports
        ensemble_report = motif_report.make_ensemble_report(cur_site_insts)
        if not reports:
            # If no results found, go back to the search page.
            message = "No search results that match your search criteria."
            return raise_validation_error(message)

        return render_to_response(
            'search_results.html',
            {'title': "",
             'description': """
             Search results can be seen as individual reports (one report per
             TF/species) or as ensemble reports (multiple TF/species). """,
             'all_cur_site_insts': [
                 pk for report in reports
                 for pk in report.get_all_cur_site_insts_ids()],
             'reports': [report.generate_browse_result_dict()
                         for report in reports],
            },
            context_instance=RequestContext(request))
    except RuntimeError as e:
        message = """
        Please select at least one TF, species and experimentaltechnique to
        search database."""
        return raise_validation_error(message)

