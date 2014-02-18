"""This file contains view functions for browsing by experimental techniques"""

from django.shortcuts import render
from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.template import RequestContext
import models
import motif_report

def browse_tech(request):
    """View function for browse by techniques. Returns all techniques, grouped
    by categories."""
    all_categories = models.ExperimentalTechniqueCategory.objects.all().order_by('name')
    binding_techniques = {}
    expression_techniques = {}
    for category in all_categories:
        # Find all techniques that belong to that category.
        # Category and techniques have n:n relationship.
        techs = models.ExperimentalTechnique.objects.filter(categories=category)
        binding_techniques[category.category_id] = techs.filter(preset_function='binding').order_by('name')
        expression_techniques[category.category_id] = techs.filter(preset_function='expression').order_by('name')

    # remove empty keys from dict
    binding_techniques = dict((x,y) for (x,y) in binding_techniques.items() if y)
    expression_techniques = dict((x,y) for (x,y) in expression_techniques.items() if y)
    all_categories = dict((x.category_id, x) for x in all_categories)
    return render_to_response('browse_tech.html',
                              {'binding_techniques': binding_techniques,
                               'expression_techniques': expression_techniques,
                               'categories': all_categories},
                              context_instance=RequestContext(request))

def get_results_tech(request, type_, id):
    """Given a technique category and an id describing that category object,
    retrieve all curation-site-instance objects that have a technique with the
    specified category."""

    techniques = None
    if type_ in ['binding', 'expression']:
        techniques = models.ExperimentalTechnique.objects.filter(preset_function=type_)
        if type_ == 'binding':
            title = 'Detection of binding'
            desc = ''
        else:
            title = 'Assessment of expression'
            desc = ''
            
    elif type_ == 'binding_category':
        category = models.ExperimentalTechniqueCategory.objects.get(category_id=id)
        techniques = models.ExperimentalTechnique.objects.filter(categories=category, preset_function='binding')
        title = category.name
        desc = category.description
    elif type_ == 'expression_category':
        category = models.ExperimentalTechniqueCategory.objects.get(category_id=id)
        techniques = models.ExperimentalTechnique.objects.filter(categories=category, preset_function='expression')
        title = category.name
        desc = category.description
    elif type_ == 'technique':
        techniques = models.ExperimentalTechnique.objects.filter(technique_id=id)
        # make sure the technique id is valid
        assert techniques.count() > 0
        title = techniques.all()[:1].get().name
        desc = techniques.all()[:1].get().description

    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        experimental_techniques=techniques
    )
    # generate all reports
    reports = motif_report.make_reports(cur_site_insts)

    return render_to_response("browse_results.html",
                              {'title': title,
                               'description': desc,
                               'reports': [report.generate_browse_result_dict() for report in reports]},
                              context_instance=RequestContext(request))
