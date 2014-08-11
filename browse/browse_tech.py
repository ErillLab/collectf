"""This file contains view functions for browsing by experimental techniques"""

from django.shortcuts import render_to_response
from django.template import RequestContext
import browse.models as models
import browse.motif_report as motif_report
from django.shortcuts import get_object_or_404, get_list_or_404

def browse_tech(request):
    """View function for browse by techniques. Returns all techniques, grouped
    by categories."""
    all_categories = models.ExperimentalTechniqueCategory.objects.all().\
                     order_by('name')
    binding_techniques = {}
    expression_techniques = {}
    for category in all_categories:
        # Find all techniques that belong to that category.
        # Category and techniques have n:n relationship.
        techs = models.ExperimentalTechnique.objects.filter(categories=category)
        binding_techniques[category.category_id] = \
                    techs.filter(preset_function='binding').order_by('name')
        expression_techniques[category.category_id] = \
                    techs.filter(preset_function='expression').order_by('name')

    # remove empty keys from dict
    binding_techniques = dict((x, y)
                              for (x, y) in binding_techniques.items() if y)
    expression_techniques = dict((x, y)
                                 for (x, y) in expression_techniques.items()
                                 if y)
    all_categories = dict((x.category_id, x) for x in all_categories)
    return render_to_response('browse_tech.html',
                              {'binding_techniques': binding_techniques,
                               'expression_techniques': expression_techniques,
                               'categories': all_categories},
                              context_instance=RequestContext(request))

def get_techniques(tech_type, tech_id):
    """Given a technique or category, return all techniques that fits under that
    category"""

    tech_objs = models.ExperimentalTechnique.objects

    name, desc, techniques = (None, None, None)
    if tech_type == 'all':
        techniques = tech_objs.all()
    elif tech_type in ['binding', 'expression']:
        techniques = tech_objs.filter(preset_function=tech_type)
        if tech_type == 'binding':
            title = 'Detection of binding'
        else:
            title = 'Assessment of expression'
        desc = ''
    elif tech_type == 'binding_category':
        category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                     category_id=tech_id)
        techniques = tech_objs.filter(categories=category,
                                      preset_function='binding')
        title = category.name
        desc = category.description
    elif tech_type == 'expression_category':
        category = get_object_or_404(models.ExperimentalTechniqueCategory,
                                     category_id=tech_id)
        techniques = tech_objs.filter(categories=category,
                                      preset_function='expression')
        title = category.name
        desc = category.description
    elif tech_type == 'technique':
        techniques = get_list_or_404(models.ExperimentalTechnique,
                                     technique_id=tech_id)
        title = techniques[0].name
        desc = techniques[0].description

    return name, desc, techniques

def get_results_tech(request, type_, id_):
    """Given a technique category and an id describing that category object,
    retrieve all curation-site-instance objects that have a technique with the
    specified category."""

    title, desc, techniques = get_techniques(type_, id_)

    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        experimental_techniques__in=techniques
    )
    # generate all reports
    reports = motif_report.make_reports(cur_site_insts)

    return render_to_response("browse_results.html",
                              dict(title=title,
                                   description=desc,
                                   reports=[report.generate_browse_result_dict()
                                            for report in reports],
                                   tax_param_type='all',
                                   tax_param=-1,
                                   tf_param_type='all',
                                   tf_param=-1,
                                   tech_param_type=type_,
                                   tech_param=id_),
                              context_instance=RequestContext(request))
