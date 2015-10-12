"""
The view functions for browsing by experimental techniques.
"""

from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404
from django.shortcuts import get_list_or_404
from django.shortcuts import render_to_response
from django.template import RequestContext

import browse.models as models
import browse.motif_report as motif_report
from browse.static_reports import get_static_reports
from browse.view_reports import view_reports_by_all_techniques
from browse.view_reports import view_reports_by_technique_category
from browse.view_reports import view_reports_by_technique

def browse_tech(request):
    """Returns all techniques, grouped by categories."""
    # Group all binding and expression techniques by their category.
    categories = models.ExperimentalTechniqueCategory.objects.all()
    binding_techs = {}
    expression_techs = {}
    for category in categories:
        techs = models.ExperimentalTechnique.objects.filter(categories=category)
        binding_techs[category.category_id] = techs.filter(
            preset_function='binding')
        expression_techs[category.category_id] = techs.filter(
            preset_function='expression')

    # remove empty keys from dict
    binding_techs = {x: y for x, y in binding_techs.items() if y}
    expression_techs = {x: y for x, y in expression_techs.items() if y}
    categories = {x.category_id: x for x in categories}
    return render_to_response('browse_tech.html',
                              {'binding_techniques': binding_techs,
                               'expression_techniques': expression_techs,
                               'categories': categories},
                              context_instance=RequestContext(request))

def get_results_all(request, function):
    """Returns motif reports of sites validated by binding or expression."""
    assert function in ['binding', 'expression']
    reports, _ = get_static_reports('experimental_technique_all_%s' % function)
    title_lookup = {'binding': 'Detection of binding',
                    'expression': 'Assessment of expression'}
    return render_to_response(
        'browse_results.html',
        {'title': title_lookup[function],
         'description': '',
         'reports': [r.generate_browse_result_dict() for r in reports],
         'combined_report_url': reverse(view_reports_by_all_techniques,
                                        args=(function,))},
        context_instance=RequestContext(request))


def get_results_category(request, category_function, object_id):
    """Returns motif reports by binding/expression category."""
    assert category_function in ['binding', 'expression']
    category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                 category_id=object_id)
    reports, _ = get_static_reports(
        'experimental_technique_category_%s' % object_id)
    return render_to_response(
        'browse_results.html',
        {'title': category.name,
         'description': category.description,
         'reports': [r.generate_browse_result_dict() for r in reports],
         'combined_report_url': reverse(
             view_reports_by_technique_category,
             args=(category_function, object_id,))},
        context_instance=RequestContext(request))

def get_results_technique(request, object_id):
    """Returns motif reports by experimental technique ID."""
    technique = get_object_or_404(models.ExperimentalTechnique,
                                  technique_id=object_id)
    reports, _ = get_static_reports('experimental_technique_%s' % object_id)
    return render_to_response(
        'browse_results.html',
        {'title': technique.name,
         'description': technique.description,
         'reports': [r.generate_browse_result_dict() for r in reports],
         'combined_report_url': reverse(view_reports_by_technique,
                                        args=(technique.technique_id,))},
        context_instance=RequestContext(request))
