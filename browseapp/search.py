from browse_base import *
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages
from browse_TF_and_species import browse_TF_and_species_selected
from browse_TF_and_species import browse_TF_and_species_selected_non_motif

def search(request):
    """Handler for search by TF/species."""
    return search_get(request) if not request.POST else search_post(request)

def search_get(request):
    binding, expression, insilico = get_techniques()
    template = {
        'TF_families': models.TFFamily.objects.all().order_by('name'),
        'phyla': models.Taxonomy.objects.filter(rank='phylum').order_by('name'),
        'binding_techniques': binding,
        'expression_techniques': expression,
        'insilico_techniques': insilico,
    }
        
    return render(request,
                  "search.html",
                  template,
                  context_instance=RequestContext(request))

def search_post(request):
    TF_input = request.POST.getlist('tf_input')
    species_input = request.POST.getlist('species_input')
    cat_input_1 = request.POST.getlist('cat_input_1')
    cat_input_2 = request.POST.getlist('cat_input_2')
    cat_input_3 = request.POST.getlist('cat_input_3')
    boolean1 = request.POST['boolean1']
    boolean2 = request.POST['boolean2']
    main_category_1 = request.POST['main_category1']
    main_category_2 = request.POST['main_category2']
    main_category_3 = request.POST['main_category3']
    
    # remove checkboxes of intermediate nodes
    TF_input = [x for x in TF_input if x!='on']
    species_input = [x for x in species_input if x!= 'on']
    cat_input1 = [x for x in cat_input_1 if x!='on']
    cat_input2 = [x for x in cat_input_2 if x!='on']
    cat_input3 = [x for x in cat_input_3 if x!='on']
 
    if not (TF_input and species_input and (cat_input_1 or cat_input2 or cat_input3)):
        message = "Please select at least one TF, species and experimental technique to search database."
        messages.add_message(request, messages.ERROR, message)
        return HttpResponseRedirect(reverse(search))
    
    TFs = models.TF.objects.filter(TF_id__in=TF_input)
    species = models.Taxonomy.objects.filter(pk__in=species_input)
    techniques1 = models.ExperimentalTechnique.objects.filter(technique_id__in=cat_input1)
    techniques2 = models.ExperimentalTechnique.objects.filter(technique_id__in=cat_input2)
    techniques3 = models.ExperimentalTechnique.objects.filter(technique_id__in=cat_input3)

    q1 = technique_list_to_Q(techniques1)
    q2 = technique_list_to_Q(techniques2)
    q3 = technique_list_to_Q(techniques3)
    # get all curation_site_instance objects satisfying those conditions
    all_csis = models.Curation_SiteInstance.objects.filter(curation__TF__in=TFs,
                                                       site_instance__genome__taxonomy__in=species)

    csis = all_csis.filter(is_motif_associated=True)
    non_motif_csis = all_csis.filter(is_motif_associated=False)

    print 'booleans', boolean1, boolean2
    if boolean1 == 'and' and boolean2 == 'and':
        csis = csis.filter(q1, q2, q3)
    elif boolean1 == 'and' and boolean2 == 'or':
        # (A and B) or C <-> (A or C) and (B or C)
        csis = csis.filter(q1|q2, q2|q3)
    elif boolean1 == 'or' and boolean2 == 'and':
        csis = csis.filter(q1|q2).filter(q3)
    elif boolean1 == 'or' and boolean2 == 'or':
        csis = csis.filter(q1|q2|q3)
    else:
        assert False, 'shouldnt be here, unhandled case'

    return render_search_results(request, csis, non_motif_csis)

def render_search_results(request, csis, non_motif_csis):
    template = search_results(csis, non_motif_csis)
    return render(request, "search_results.html", template, context_instance=RequestContext(request))

def search_results(csis, non_motif_csis):
    values = csis.values('curation__TF', 'site_instance__genome__taxonomy')\
             .distinct()\
             .order_by('curation__TF__name', 'site_instance__genome__taxonomy__name')
    
    all_reports = []
    for val in values:
        TF_id = val['curation__TF']
        species_id = val['site_instance__genome__taxonomy']
        TF = models.TF.objects.get(TF_id=TF_id)
        species = models.Taxonomy.objects.get(pk=species_id)
        filtered_csis = csis.filter(curation__TF=TF, site_instance__genome__taxonomy=species)
        # use non-motif sites only if there is motif-associated data
        filtered_non_motif_csis = non_motif_csis.filter(curation__TF=TF, site_instance__genome__taxonomy=species)
        if filtered_csis:
            all_reports.append({
                'TF_name': TF.name,
                'species_name': species.name,
                'csi_list': [csi.pk for csi in filtered_csis.iterator()],
                'non_motif_csis': [ncsi.pk for ncsi in filtered_non_motif_csis.iterator()]
            })

    return {'reports': all_reports}
    
def technique_list_to_Q(techniques):
    # prepare filters
    q = Q(curation__curation_id=-9999)
    for t in techniques:
        q = q | Q(curation__experimental_techniques=t)
    return q

def get_techniques():
    binding_techniques = {}
    expression_techniques = {}
    insilico_techniques = {}
    all_categories = models.ExperimentalTechniqueCategory.objects.all()
    for category in all_categories:
        techs = models.ExperimentalTechnique.objects.filter(categories=category).order_by('name')
        binding_techniques[category.name] = techs.filter(preset_function='binding')
        expression_techniques[category.name] = techs.filter(preset_function='expression')
        insilico_techniques[category.name] = techs.filter(preset_function='insilico')
    # remove empty keys
    binding_techniques = dict((x,y) for (x,y) in binding_techniques.items() if y)
    expression_techniques = dict((x,y) for (x,y) in expression_techniques.items() if y)
    insilico_techniques = dict((x,y) for (x,y) in insilico_techniques.items() if y)

    return binding_techniques, expression_techniques, insilico_techniques



