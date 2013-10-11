from browse_base import *
import baseapp.templatetags.TF_desc_htmlify as htmlify
from django.forms.models import model_to_dict
from django.core.urlresolvers import reverse
from browse_TF_and_species import browse_TF_and_species_selected
import Queue
from search import group_search_results

def browse_TF(request):
    TF_families = models.TFFamily.objects.all().order_by('name')
    return render_to_response('browse_TF.html',
                              {'TF_families': TF_families},
                              context_instance=RequestContext(request))

def browse_TF_all_reports_ajax(request, t, id):
    # Return links to all reports
    # type can be TF or TF_family
    if t=='TF':
        TFs = [models.TF.objects.get(TF_id=id)]
        TF_name = TFs[0].name
        desc = TFs[0].description
    elif t=='TF_family':
        TF_family = models.TFFamily.objects.get(TF_family_id=id)
        TFs = models.TF.objects.filter(family=TF_family).all()
        desc = TF_family.description
        TF_name = TF_family.name
        
    motif_csis = models.Curation_SiteInstance.objects.filter(curation__TF__in=TFs,
                                                             is_motif_associated=True)
    non_motif_csis = models.Curation_SiteInstance.objects.filter(curation__TF__in=TFs,
                                                                 is_motif_associated=False)
    all_reports = group_search_results(motif_csis, non_motif_csis)['reports']

    assert TF_name
    return render_to_response("browse_tab.html",
                              {'title': TF_name,
                               'description': desc,
                               'reports': all_reports,
                               'view_all_csis': [csi_pk for report in all_reports for csi_pk in report['csi_list']],
                               'view_all_ncsis': [ncsi_pk for report in all_reports for ncsi_pk in report['non_motif_csis']]},
                              context_instance=RequestContext(request))
    
def browse_tax(request):
    all_species_ids = models.Curation_SiteInstance.objects.values_list('site_instance__genome__taxonomy__pk', flat=True).distinct()
    all_species = models.Taxonomy.objects
    # get all phylums
    taxonomy = {'phyla': all_species.filter(rank='phylum')}
    return render_to_response('browse_SP.html',
                              {'taxonomy': taxonomy},
                              context_instance=RequestContext(request))
    
def browse_tax_all_reports_ajax(request, id):
    all_sp = []
    Q = Queue.Queue()
    Q.put(models.Taxonomy.objects.get(pk=id))
    while not Q.empty():
        s = Q.get()
        children =  s.taxonomy_set.all()
        if children:
            for c in children:
                Q.put(c)
        else:
            all_sp.append(s)

    motif_csis = models.Curation_SiteInstance.objects.filter(site_instance__genome__taxonomy__in=all_sp, is_motif_associated=True)
    non_motif_csis = models.Curation_SiteInstance.objects.filter(site_instance__genome__taxonomy__in=all_sp, is_motif_associated=False)
    all_reports = group_search_results(motif_csis, non_motif_csis)['reports']
    return render_to_response("browse_tab.html",
                              {'title': models.Taxonomy.objects.get(pk=id).name,
                               'description': '',
                               'reports': all_reports},
                              context_instance=RequestContext(request))

def browse_techniques(request):
    all_categories = models.ExperimentalTechniqueCategory.objects.all().order_by('name')
    binding_techniques = {}
    expression_techniques = {}
    for category in all_categories:
        # find all techniques that belong to that category
        # category and techniques have n:n relationship
        techs = models.ExperimentalTechnique.objects.filter(categories=category)
        binding_techniques[category.category_id] = techs.filter(preset_function='binding').order_by('name')
        expression_techniques[category.category_id] = techs.filter(preset_function='expression').order_by('name')

    # remove empty keys from dict
    binding_techniques = dict((x,y) for (x,y) in binding_techniques.items() if y)
    expression_techniques = dict((x,y) for (x,y) in expression_techniques.items() if y)
    all_categories = dict((x.category_id, x) for x in all_categories)
    return render_to_response('browse_TECH.html',
                              {'binding_techniques': binding_techniques,
                               'expression_techniques': expression_techniques,
                               'categories': all_categories},
                              context_instance=RequestContext(request))

def browse_techniques_all_reports_ajax(request, type, id):
    techniques = None
    if type in ['binding', 'expression']:
        techniques = models.ExperimentalTechnique.objects.filter(preset_function=type)
        if type=='binding':
            title = 'Detection of binding'
            desc = ''
        else:
            title = 'Assessment of expression'
            desc = ''
            
    elif type=='binding_category':
        category = models.ExperimentalTechniqueCategory.objects.get(category_id=id)
        techniques = models.ExperimentalTechnique.objects.filter(categories=category, preset_function='binding')
        title = category.name
        desc = category.description
    elif type=='expression_category':
        category = models.ExperimentalTechniqueCategory.objects.get(category_id=id)
        techniques = models.ExperimentalTechnique.objects.filter(categories=category, preset_function='expression')
        title = category.name
        desc = category.description
    elif type=='technique':
        techniques = models.ExperimentalTechnique.objects.filter(technique_id=id)
        assert techniques.count() > 0
        title = techniques.all()[:1].get().name
        desc = techniques.all()[:1].get().description

    motif_csis = models.Curation_SiteInstance.objects.filter(curation__experimental_techniques__in=techniques, is_motif_associated=True)
    non_motif_csis = models.Curation_SiteInstance.objects.filter(curation__experimental_techniques__in=techniques, is_motif_associated=False)
    all_reports = group_search_results(motif_csis, non_motif_csis)['reports']
    return render_to_response("browse_tab.html",
                              {'title': title,
                               'description': desc,
                               'reports': all_reports},
                              context_instance=RequestContext(request))
