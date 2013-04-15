from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.shortcuts import render
from django.template import RequestContext
from collectfapp.bioutils import weblogo
from base64 import b64encode
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db.models import Q

import models
import fetch
import forms

strain_orders = {}  # add a new field for strain table in database which holds order info

@login_required
def view_all_curations(request):
    """Handler function to see all curations at once"""
    template_vals = {"curations": fetch.get_all_curations()}
    return render_to_response("curation_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))

@login_required
def view_all_publications(request):
    """Handler function to see all publications in the database"""
    template_vals = {"publications": fetch.get_all_publications()}
    return render_to_response("publication_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))

def browse_get(request):
    """Return empty browse form"""
    form = forms.BrowseForm()
    return render(request,
                  "browse.html",
                  {'form': forms.BrowseForm(),},
                  context_instance=RequestContext(request))

def browse_post(request):
    """Process form for browsing"""
    form = forms.BrowseForm(request.POST)
    if form.is_valid():
        TF = form.cleaned_data['TF']
        species = form.cleaned_data['species']
        experimental_techniques = form.cleaned_data['techniques']
    return get_sites_by_TF_species(request, TF, species, experimental_techniques)

def browse_post_TF_sp(request, TF_id, species_id):
    """Handle Http requests with TF_id and species_id"""
    TF = fetch.get_TF_by_id(TF_id)
    species = fetch.get_species_by_id(species_id)
    return get_sites_by_TF_species(request, TF, species,
                                   experimental_techniques=models.ExperimentalTechnique.objects.all())
    
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

def get_sites_by_TF_species(request, TF, species, experimental_techniques):
    """Given a TF, species and a list of experimental techniques, query the database
    to search sites for the particular TF and species, filtering out not selected
    experimental techniques"""
    # see collectfapp.models for description of model objects
    # get all Curation_SiteInstance objects
    curation_site_instances = fetch.get_curation_site_instances(TF, species)
    print len(curation_site_instances)
    # filter curation_site_instances based on used experimental techniques
    if experimental_techniques:
        x = Q()
        for exp_technique in experimental_techniques:
            x = x | Q(curation__experimental_techniques=exp_technique)
            
        curation_site_instances = curation_site_instances.filter(x).distinct()
    else:
        curation_site_instances = []

    print len(curation_site_instances)
    # group them by site instance?
    site_curation_dict, site_regulation_dict = group_curation_site_instances(curation_site_instances)

    site_sequences = set(csi.site_instance for csi in curation_site_instances)

    # create weblogo for the list of sites
    #weblogo_data = weblogo_uri(map(lambda x: x.seq, site_sequences))
    # if there is no site, message
    if not site_curation_dict and not site_regulation_dict:
        messages.info(request, "No site found for transcription factor %s in the genome of %s." % (TF.name,
                                                                                                   species.name))

    response_dict = {'form': forms.BrowseForm(initial={'TF': TF,
                                                       'species': species,
                                                       'experimental_techniques': experimental_techniques}),
                     'TF': TF,
                     'species': species,
                     'site_curation_dict': site_curation_dict,
                     'site_regulation_dict': site_regulation_dict,
                     #'weblogo_image_data': weblogo_data
                     }
    return render(request, "browse.html", response_dict, context_instance=RequestContext(request))
    
def browse(request):
    """Handler function to browse database"""
    if not request.POST:
        return browse_get(request)

    return browse_post(request)

def weblogo_uri(sequences):
    """Generate the weblogo and make it ready for direct embed into response html"""
    image_data = weblogo(sequences)
    encoded = b64encode(image_data)
    mime = "image/png"
    return ("data:" + mime + ';' + "base64," + encoded)

def export_sites(request):
    """Given a list of sites, report FASTA/CSV file containing sites for particular TF
    and species"""
    if 'fasta' in request.POST: export_format = 'fasta'
    elif 'csv'in request.POST: export_format = 'csv'
    
    assert export_format in ['fasta', 'csv']
    site_ids = request.POST.getlist('site_id')
    sites = models.SiteInstance.objects.filter(site_id__in=site_ids)
    filename = "sites.fasta" if export_format == 'fasta' else "sites.csv"
    # set HttpResponse stuff
    response = HttpResponse(content_type='application/download')
    response['Content-Disposition'] = 'attachment;filename="%s"' % filename
    # write all sites to file
    for site in sites:
        response.write(site.to_fasta() if export_format=='fasta' else site.to_csv())
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
            num_site_instances = curation_site_instances.values_list('site_instance', flat=True).distinct().count()
            num_curations = curation_site_instances.values_list('curation', flat=True).distinct().count()
            curation_stats[TF.name][sp.name] = num_curations
            curation_stats2[TF.name][sp.name]= num_site_instances

    response_dict = dict(curation_stats=curation_stats, curation_stats2=curation_stats2,
                         TFs = [tf.name for tf in all_TFs], species = [sp.name for sp in all_species])
    return render(request, "database_stats.html", response_dict, context_instance=RequestContext(request))
    

def browse_by_TF_main(request):
    """Handler for browse by TF request"""
    TFs = fetch.get_all_TFs()
    TF_families = fetch.get_all_TF_families()
    response_dict = {'TF_families': TF_families}
    return render(request, "browse_tf_main.html", response_dict, context_instance=RequestContext(request))

def browse_by_TF_family(request, TF_family_id):
    """Handler for browse TF family"""
    TF_family = fetch.get_TF_family_by_id(TF_family_id)
    TFs = fetch.get_TFs_by_family(TF_family)
    response_dict = {'TF_family': TF_family, 'TFs': TFs}
    return render(request, "browse_tf_family.html", response_dict, context_instance=RequestContext(request))


def browse_by_TF(request, TF_id):
    """Handler for browse TF"""
    TF = fetch.get_TF_by_id(TF_id)
    # fetch species that have curation data on this TF
    curation_site_instances = models.Curation_SiteInstance.objects.filter(curation__TF=TF)
    species_ids = curation_site_instances.values_list('site_instance__genome__strain', flat=True).distinct()
    species = models.Strain.objects.filter(taxonomy_id__in=species_ids)

    num_site_instances = {} # dictionary of num_site_instances by species_id
    num_curations = {}      # dictionary of num_curations by species_id
    # get number of sites and curations for each species
    for sp in species:
        # get filtered curation_site_instances
        filtered_csi = curation_site_instances.filter(site_instance__genome__strain=sp)
        num_site_instances[sp.taxonomy_id] = filtered_csi.values_list('site_instance', flat=True).distinct().count()
        num_curations[sp.taxonomy_id] = filtered_csi.values_list('curation', flat=True).distinct().count()

    response_dict = {'TF': TF,
                     'species': species,
                     'num_site_instances': num_site_instances,
                     'num_curations': num_curations}
    return render(request, "browse_tf.html", response_dict, context_instance=RequestContext(request))
    

