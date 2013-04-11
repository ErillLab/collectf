from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.shortcuts import render
from django.template import RequestContext
from collectfapp.bioutils import weblogo
from base64 import b64encode

import models
import fetch
import forms

strain_orders = {}  # add a new field for strain table in database which holds order info

def view_all_curations(request):
    """Handler function to see all curations at once"""
    template_vals = {
        "curations": fetch.get_all_curations()
    }
    return render_to_response("curation_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))

def view_all_publications(request):
    """Handler function to see all publications in the database"""
    template_vals = {
        "publications": fetch.get_all_publications()
    }
    return render_to_response("publication_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))

def browse_get(request):
    form = forms.BrowseForm()
    return render(request,
                  "browse.html",
                  {'form': forms.BrowseForm(),},
                  context_instance=RequestContext(request))

def browse_post(request):
    form = forms.BrowseForm(request.POST)
    if form.is_valid():
        TF_id = form.cleaned_data['TF']
        species_id = form.cleaned_data['species']
        experimental_techniques = form.cleaned_data['techniques']

    return browse_post_helper(request, TF_id, species_id, experimental_techniques)

def browse_post_TF_sp(request, TF_id, species_id):
    return browse_post_helper(request, TF_id, species_id, experimental_techniques=None)
    

def browse_post_helper(request, TF_id, species_id, experimental_techniques):

    # see collectfapp.models for description of model objects
    # get all Curation_SiteInstance objects
    curation_site_instances = fetch.get_curation_site_instances(TF_id, species_id)
    # filter curation_site_instances based on used experimental techniques
    # filter goes here
    
    # group them by site instance?
    site_curation_dict, site_regulation_dict = group_curation_site_instances(curation_site_instances)

    site_sequences = set(csi.site_instance for csi in curation_site_instances)
    weblogo_data = weblogo_uri(map(lambda x: x.seq, site_sequences))

    response_dict = {
        'form': forms.BrowseForm(initial={'TF': TF_id,
                                          'species': species_id,
                                          'experimental_techniques': experimental_techniques}),
        'site_curation_dict': site_curation_dict,
        'site_regulation_dict': site_regulation_dict,
        'weblogo_image_data': weblogo_data,
        }
    return render(request, "browse.html", response_dict, context_instance=RequestContext(request))
    
def browse(request):
    """Handler function to browse database"""
    if not request.POST:
        return browse_get(request)

    return browse_post(request)

def display_weblogo(request, image_data):
    return HttpResponse(image_data, mimetype="image/png")

def weblogo_uri(sequences):
    """Generate the weblogo and make it ready for direct embed into response html"""
    image_data = weblogo(sequences)
    encoded = b64encode(image_data)
    mime = "image/png"
    return ("data:" + mime + ';' + "base64," + encoded)

def report_FASTA():
    """Given a list of sites, report FASTA file containing sites for particular TF
    and species"""
    # get TF and species from id
    TF = fetch.get_TF_by_id(TF_id)
    species = fetch.get_species_by_id(species_id) 

    # set HttpResponse stuff
    response = HttpResponse(content_type='application/download')
    response['Content-Disposition'] = 'attachment;filename="%s_%s.fasta"' % (TF.name, species.name)
    
    # write all sites to file
    for curation_site_instance in fetch.get_curations(TF_id, species_id):
        response.write(curation_site_instance.site_instance.to_fasta())
        
    return response
    
def curation_stats(request):
    all_TFs = fetch.get_all_TFs()
    all_species = fetch.get_all_species()
    # number of curations by TF and site
    curation_stats = {}  # dictionary of number of curations by TF and species
    curation_stats2 = {} # dictionary of number of sites by TF and species
    for TF in all_TFs:
        curation_stats[TF.name] = {} # dict[species]
        curation_stats2[TF.name] = {}
        for sp in all_species:
            curation_site_instances = fetch.get_curation_site_instances(TF, sp)
            curations = set(csi.curation for csi in curation_site_instances)
            sites = set(csi.site_instance for csi in curation_site_instances)
            curation_stats[TF.name][sp.name] = len(curations)  # num of curations
            curation_stats2[TF.name][sp.name]= len(sites)      # num of sites


    return render(request,
                  "database_stats.html",
                  dict(
                      curation_stats = curation_stats,
                      curation_stats2 = curation_stats2,
                      TFs = [tf.name for tf in all_TFs],
                      species = [sp.name for sp in all_species]
                      ),
                  context_instance=RequestContext(request))
    


def browse_by_species_main(request):
    """Handler for browse by species request"""
    species = fetch.get_all_species()
    for sp in species:
        if sp.name not in strain_orders:
            print '********************************************************************************'
            print 'entrez taxonomy', sp.name, sp.taxonomy_id
            order_name = fetch.get_lineage(sp.taxonomy_id)
            strain_orders[sp.name] = order_name

    orders = {}
    for sp in species:
        orders[strain_orders[sp.name]] = sp
            
    response_dict = {'orders': orders}
    return render(request,
                  "browse_sp_main.html",
                  response_dict,
                  context_instance=RequestContext(request))


def browse_by_TF_main(request):
    """Handler for browse by TF request"""
    TFs = fetch.get_all_TFs()
    TF_families = fetch.get_all_TF_families()
    # group TFs by TF family    
    response_dict = {'TF_families': TF_families}
    return render(request,
                  "browse_tf_main.html",
                  response_dict,
                  context_instance=RequestContext(request))

def browse_by_TF_family(request, TF_family_id):
    """Handler for browse TF family"""
    TF_family = fetch.get_TF_family_by_id(TF_family_id)
    TFs = fetch.get_TFs_by_family(TF_family)
    response_dict = {'TF_family': TF_family,
                     'TFs': TFs}
    return render(request, "browse_tf_family.html", response_dict, context_instance=RequestContext(request))


def browse_by_TF(request, TF_id):
    """Handler for browse TF"""
    TF = fetch.get_TF_by_id(TF_id)
    # fetch species that have curation data on this TF
    species_ids = models.Curation_SiteInstance.objects.filter(curation__TF=TF).values_list('site_instance__genome__strain')
    species_ids = list(set(species_ids))
    species = models.Strain.objects.filter(taxonomy_id__in=map(lambda x:x[0], species_ids))
    response_dict = {'TF': TF,
                     'species': species}
    return render(request, "browse_tf.html", response_dict, context_instance=RequestContext(request))
    

def group_curation_site_instances(curation_site_instances):
    """Group curation_site_instance objects by site_instance"""
    site_curation_dict = dict((csi.site_instance,[]) for csi in curation_site_instances)
    site_regulation_dict = dict((csi.site_instance,[]) for csi in curation_site_instances)
    for csi in curation_site_instances:
        s = csi.site_instance
        site_curation_dict[s].append(csi.curation)  # insert curation
        for reg in csi.regulation_set.all():
            # check if regulated gene is already in the list (from another curation)
            same_gene_reg = [r for r in site_regulation_dict[s] if r.gene == reg.gene]
            if not same_gene_reg:
                site_regulation_dict[s].append(reg)
            elif same_gene_reg[0].evidence_type == "inferred" and reg.evidence_type == "exp_verified":
                same_gene_reg[0].evidence_type = "exp_verified"
        
    return site_curation_dict, site_regulation_dict
