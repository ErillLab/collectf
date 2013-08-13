from browse_base import *
import baseapp.templatetags.TF_desc_htmlify as htmlify
from django.forms.models import model_to_dict
from django.core.urlresolvers import reverse
from browse_TF_and_species import browse_TF_and_species_selected
import Queue

def browse_TF(request):
    TF_families = models.TFFamily.objects.all().order_by('name')
    return render_to_response('browse_TF.html',
                              {'TF_families': TF_families},
                              context_instance=RequestContext(request))

def browse_TF_all_reports_json(request, t, id):
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
        
    TF_sp = models.Curation_SiteInstance.objects.filter(curation__TF__in=TFs, is_motif_associated=True)\
            .order_by('curation__TF__name', 'site_instance__genome__taxonomy__name')\
            .values('curation__TF__name',
                   'site_instance__genome__taxonomy__name',
                   'curation__TF__TF_id',
                   'site_instance__genome__taxonomy__pk').distinct()

    all_reports = [{
        'TF_name': elm['curation__TF__name'],
        'species_name': elm['site_instance__genome__taxonomy__name'],
        'TF_id': elm['curation__TF__TF_id'],
        'species_id': elm['site_instance__genome__taxonomy__pk'],
        'report_link': reverse(browse_TF_and_species_selected, kwargs={'TF_param': 'TF',
                                                                       'species_param': 'species',
                                                                       'TF_ids': elm['curation__TF__TF_id'],
                                                                       'species_ids': elm['site_instance__genome__taxonomy__pk']}),
        } for elm in TF_sp]

    TF_ids = list(set([str(r['TF_id']) for r in all_reports]))
    species_ids = list(set([str(r['species_id']) for r in all_reports]))
    view_all_report_link = None
    if TF_ids and species_ids:
        view_all_report_link = reverse(browse_TF_and_species_selected, kwargs={'TF_param': 'TF',
                                                                               'species_param': 'species',
                                                                               'TF_ids': ','.join(TF_ids),
                                                                               'species_ids': ','.join(species_ids)})
    return HttpResponse(json.dumps({'description': htmlify.htmlify(desc),
                                    'TF_name': TF_name,
                                    'list': all_reports,
                                    'view_all_report_link': view_all_report_link}),
                        mimetype="application/json")

    
def browse_tax(request):
    all_species_ids = models.Curation_SiteInstance.objects.values_list('site_instance__genome__taxonomy__pk', flat=True).distinct()
    all_species = models.Taxonomy.objects
    # get all phylums
    """
    roots = set()
    for sp in all_species:
        while sp.parent:
            sp = sp.parent
        roots.add(sp)

    roots = sorted(list(roots), key=lambda x: x.name)
    """
    taxonomy = {
        'phyla': all_species.filter(rank='phylum')
    }
    return render_to_response('browse_SP.html',
                              {'taxonomy': taxonomy},
                              context_instance=RequestContext(request))
    
def browse_tax_all_reports_json(request, id):
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
    
    csis = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in=all_sp,
        is_motif_associated=True)\
        .order_by('curation__TF__name', 'site_instance__genome__taxonomy__name')\
        .values('curation__TF__TF_id',
                'curation__TF__name',
                'curation__TF_instance__protein_accession',
                'site_instance__genome__taxonomy__pk',
                'site_instance__genome__taxonomy__name').distinct()
    
    all_reports = [{
        'TF_name': elm['curation__TF__name'],
        'TF_accession': elm['curation__TF_instance__protein_accession'],
        'species_name': elm['site_instance__genome__taxonomy__name'],
        'TF_id': elm['curation__TF__TF_id'],
        'species_id': elm['site_instance__genome__taxonomy__pk'],
        'report_link': reverse(browse_TF_and_species_selected, kwargs={'TF_param': 'TF',
                                                                       'species_param': 'species',
                                                                       'TF_ids': elm['curation__TF__TF_id'],
                                                                       'species_ids': elm['site_instance__genome__taxonomy__pk']}),
        } for elm in csis]

    TF_ids = list(set([str(r['TF_id']) for r in all_reports]))
    species_ids = list(set([str(r['species_id']) for r in all_reports]))
    view_all_report_link = None
    if TF_ids and species_ids:
        view_all_report_link = reverse(browse_TF_and_species_selected, kwargs={'TF_param': 'TF',
                                                                               'species_param': 'species',
                                                                               'TF_ids': ','.join(TF_ids),
                                                                               'species_ids': ','.join(species_ids)})
    return HttpResponse(json.dumps({'description': 'desc goes here',
                                    'list': all_reports,
                                    'view_all_report_link': view_all_report_link}),
                        mimetype="application/json")

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

    print binding_techniques
    return render_to_response('browse_TECH.html',
                              {'binding_techniques': binding_techniques,
                               'expression_techniques': expression_techniques,
                               'categories': all_categories},
                              context_instance=RequestContext(request))
    

def browse_techniques_all_reports_json(request, type, id):
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

        
    csis = models.Curation_SiteInstance.objects.filter(curation__experimental_techniques__in=techniques,
                                                       is_motif_associated=True)\
           .order_by('curation__TF__name', 'site_instance__genome__taxonomy__name')\
           .values('curation__TF_instance__protein_accession',
                   'curation__TF__name',
                   'curation__TF__TF_id',
                   'site_instance__genome__taxonomy__pk',
                   'site_instance__genome__taxonomy__name').distinct()
    
    all_reports = [{
        'TF_name': elm['curation__TF__name'],
        'TF_accession': elm['curation__TF_instance__protein_accession'],
        'species_name': elm['site_instance__genome__taxonomy__name'],
        'TF_id': elm['curation__TF__TF_id'],
        'species_id': elm['site_instance__genome__taxonomy__pk'],
        'report_link': reverse(browse_TF_and_species_selected, kwargs={'TF_param': 'TF',
                                                                       'species_param': 'species',
                                                                       'TF_ids': elm['curation__TF__TF_id'],
                                                                       'species_ids': elm['site_instance__genome__taxonomy__pk']}),
        } for elm in csis]

    TF_ids = list(set([str(r['TF_id']) for r in all_reports]))
    species_ids = list(set([str(r['species_id']) for r in all_reports]))
    view_all_report_link = None
    if TF_ids and species_ids:
        view_all_report_link = reverse(browse_TF_and_species_selected, kwargs={'TF_param': 'TF',
                                                                               'species_param': 'species',
                                                                               'TF_ids': ','.join(TF_ids),
                                                                               'species_ids': ','.join(species_ids)})

    return HttpResponse(json.dumps({'description': desc,
                                    'title': title,
                                    'list': all_reports,
                                    'view_all_report_link': view_all_report_link}),
                        mimetype="application/json")
