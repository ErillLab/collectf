from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.shortcuts import render
from django.template import RequestContext
from collectfapp.bioutils import weblogo
from base64 import b64encode
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db.models import Q

import utils
import lasagna
import models
import fetch
import forms
import json


# View Functions

@login_required
def view_all_curations(request):
    """Handler function to see all curations at once. This function renders the page
    with the list of all curations in the database"""
    template_vals = {"curations": fetch.get_all_curations()}
    return render_to_response("curation_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))

@login_required
def view_all_publications(request):
    """Handler function to see all publications in the database. This is for internal
    use, to see all publications in the database."""
    template_vals = {"publications": fetch.get_all_publications()}
    return render_to_response("publication_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))

def browse(request):
    """Handler function to browse database. Either calls get or post view function."""
    if not request.POST:
        return browse_get(request)
    return browse_post(request)

def browse_get(request):
    """Handler for main browse page. Renders the page with the form which enables
    user to query the database."""
    return render(request, "browse.html",
                  make_browse_response_dict(),
                  context_instance=RequestContext(request))

def browse_post(request):
    """Process form to get selected TF, species and experimental techniques. Retrieve
    binding sites for selected TF and species from database, filter them by
    experimental techniques, return to the user."""
    TF = fetch.get_TF_by_id(request.POST['TF'])
    species = fetch.get_species_by_id(request.POST['species'])

    experimental_techniques_1 = request.POST['tech1']
    experimental_techniques_2 = request.POST['tech2']
    experimental_techniques_3 = request.POST['tech3']
    boolean1 = request.POST['boolean1']
    boolean2 = request.POST['boolean2']
    assert boolean1 in ['and', 'or'] and boolean2 in ['and', 'or']
    assert 'tech1' in request.POST and 'tech2' in request.POST and 'tech3' in request.POST

    print request.POST
    q1 = techniques_JSON_to_Q(experimental_techniques_1)
    q2 = techniques_JSON_to_Q(experimental_techniques_2)
    q3 = techniques_JSON_to_Q(experimental_techniques_3)

    #  get curation objects
    curation_site_instances = fetch.get_curation_site_instances(TF, species)
    # filter them by experiemntal techniques
    print len(curation_site_instances)
    if boolean1 == 'and' and boolean2 == 'and':
        curation_site_instances = curation_site_instances.filter(q1).filter(q2).filter(q3)
    elif boolean1 == 'and' and boolean2 == 'or':
        # (A and B) or C <-> (A or C) and (B or C)
        curation_site_instances = curation_site_instances.filter(q1 | q3).filter(q2 | q3)
    elif boolean1 == 'or' and boolean2 == 'and':
        curation_site_instances = curation_site_instances.filter(q1 | q2).filter(q3)
    elif boolean1 == 'or' and boolean2 == 'or':
        curation_site_instances = curation_site_instances.filter(q1 | q2 | q3)
    else:
        assert False, "shouldn't be here, query"
    
    print len(curation_site_instances)

    return get_sites_by_TF_species(request, TF, species, curation_site_instances)


def techniques_JSON_to_Q(JSON_string):
    """Given JSON string received from the form, parse the techniques (names), get
    the ids for them, build the Q object for filtering"""
    j = json.loads(JSON_string)

    techniques = map(lambda x: map(lambda y: y['key'], x['values']), j)
    techniques = [t for grp in techniques for t in grp] # flatten list
    print techniques
    q = Q(curation__curation_id=-9999)
    
    for t in techniques:
        technique_object = models.ExperimentalTechnique.objects.get(name=t)
        print technique_object
        q = q | Q(curation__experimental_techniques=technique_object)
    return q

def browse_post_TF_sp(request, TF_id, species_id):
    """Handle Http requests with TF_id and species_id"""
    TF = fetch.get_TF_by_id(TF_id)
    species = fetch.get_species_by_id(species_id)
    curation_site_instances = fetch.get_curation_site_instances(TF, species)

    return get_sites_by_TF_species(request, TF, species, curation_site_instances)

def browse_by_TF_main(request):
    """Handler for browse by TF request"""
    TFs = fetch.get_all_TFs()
    TF_families = fetch.get_all_TF_families()
    response_dict = {'TF_families': TF_families}
    return render(request, "browse_tf_main.html", response_dict,
                  context_instance=RequestContext(request))

def browse_by_TF_family(request, TF_family_id):
    """Handler for browse TF family"""
    TF_family = fetch.get_TF_family_by_id(TF_family_id)
    TFs = fetch.get_TFs_by_family(TF_family)
    response_dict = {'TF_family': TF_family, 'TFs': TFs}
    return render(request, "browse_tf_family.html", response_dict,
                  context_instance=RequestContext(request))

def browse_by_TF(request, TF_id):
    """Handler for browse TF"""
    TF = fetch.get_TF_by_id(TF_id)
    # fetch species that have curation data on this TF
    csi = models.Curation_SiteInstance.objects.filter(curation__TF=TF)
    species_ids = csi.values_list('site_instance__genome__strain', flat=True).distinct()
    species = models.Strain.objects.filter(taxonomy_id__in=species_ids)

    num_site_instances = {} # dictionary of num_site_instances by species_id
    num_curations = {}      # dictionary of num_curations by species_id
    # get number of sites and curations for each species
    for sp in species:
        # get filtered curation_site_instances
        filtered_csi = csi.filter(site_instance__genome__strain=sp)
        num_site_instances[sp.taxonomy_id] = filtered_csi.values_list('site_instance', flat=True).distinct().count()
        num_curations[sp.taxonomy_id] = filtered_csi.values_list('curation', flat=True).distinct().count()

    response_dict = {'TF': TF,
                     'species': species,
                     'num_site_instances': num_site_instances,
                     'num_curations': num_curations}
    return render(request, "browse_tf.html", response_dict,
                  context_instance=RequestContext(request))

def browse_by_site(request, dbxref_id):
    """Handler for browsing site instances"""
    site_instance_id = utils.dbxref2id(dbxref_id)
    site_instance = models.SiteInstance.objects.get(site_id=site_instance_id)
    response_dict = {'site_instance': site_instance}
    return render(request, "browse_site.html", response_dict,
                  context_instance=RequestContext(request))

def curation_stats(request):
    """Handler for curation statistics page. Count the number of curations/sites for
    each TF and species in the database, pass dictionary to the HTML template"""
    all_TFs = fetch.get_all_TFs()
    all_species = fetch.get_all_species()
    # number of curations by TF and site
    curation_stats = {}  # dictionary of number of curations by TF and species
    curation_stats2 = {} # dictionary of number of sites by TF and species
    for TF in all_TFs:
        curation_stats[TF.name] = {} # dict[species]
        curation_stats2[TF.name] = {} 
        for sp in all_species:
            csi = fetch.get_curation_site_instances(TF, sp) # curation_site_instances
            num_sites = csi.values_list('site_instance', flat=True).distinct().count()
            num_curations = csi.values_list('curation', flat=True).distinct().count()
            curation_stats[TF.name][sp.name] = num_curations
            curation_stats2[TF.name][sp.name]= num_sites

    response_dict = dict(
        number_TFs = len(all_TFs),
        number_species = len(all_species),
        number_curations = len(models.Curation.objects.all()),
        number_sites = len(models.SiteInstance.objects.all()),
        number_publications = len(models.Publication.objects.all()),
        pub_completed = '%.1f' % (len(models.Publication.objects.filter(curation_complete=True)) * 100.0 /
                                  len(models.Publication.objects.all())),
        curation_stats=curation_stats,
        curation_stats2=curation_stats2,
        TFs = [tf.name for tf in all_TFs],
        species = [sp.name for sp in all_species])
    
    return render(request, "database_stats.html", response_dict,
                  context_instance=RequestContext(request))

def export_sites(request):
    """Given a list of sites, report FASTA/CSV file containing sites for particular
    TF and species"""
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
    if export_format == 'csv':
        response.write('\t'.join(['genome', 'site_start', 'site_end', 'site_strand', 'site_seq', 'regulated_genes', 'references']))
        response.write('\n')
    for site in sites:
        response.write(site.to_fasta() if export_format=='fasta' else site.to_csv())
        if export_format == 'csv':
            """
            techniques = []
            for csi in site.curation_siteinstance_set.all():
                techniques.append(', '.join(t.name for t in csi.curation.experimental_techniques.all()) + '(%s)' % csi.curation.publication.pmid)
            response.write('\t' + '\t'.join(techniques))
            """
            regulations = models.Regulation.objects.filter(curation_site_instance__site_instance=site)

            response.write('\t' + ','.join(reg.gene.name for reg in regulations))
            response.write('\n')
    return response


# View helper functions

def make_browse_response_dict():
    """For several browse handlers (browse.GET, browse.POST, ..), the form to search
    database is always displayed on top. The data that go in this form is always
    same. Instead of preparing same data in all view functions, create a dict in this
    one and call this function from them."""
    # group techniques
    # binding
    techniques = [t for t in models.ExperimentalTechnique.objects.filter(category__main_category='binding')]
    binding_techniques = {} # dictionary of (subcategory, list of techniques in that subcategory)
    for t in techniques:
        cat = t.category
        binding_techniques[cat] = binding_techniques.get(cat, []) + [t]
    # induction
    techniques = [t for t in models.ExperimentalTechnique.objects.filter(category__main_category='induction')]
    induction_techniques = {}
    for t in techniques:
        cat = t.category
        induction_techniques[cat] = induction_techniques.get(cat, []) + [t]
    # insilico
    techniques = [t for t in models.ExperimentalTechnique.objects.filter(category__main_category='insilico')]
    insilico_techniques = {}
    for t in techniques:
        cat = t.category
        insilico_techniques[cat] = insilico_techniques.get(cat, []) + [t]

    return dict(TFs=models.TF.objects.all(),
                species=models.Strain.objects.all(),
                binding_techniques=binding_techniques,
                induction_techniques=induction_techniques,
                insilico_techniques=insilico_techniques)
    
   
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

def get_sites_by_TF_species(request, TF, species, curation_site_instances):    
    # group them by site instance?
    site_curation_dict, site_regulation_dict = group_curation_site_instances(curation_site_instances)
    site_sequences = set(csi.site_instance for csi in curation_site_instances)

    # if there is no site, message
    if not site_curation_dict and not site_regulation_dict:
        msg = "No site found for transcription factor %s in the genome of %s." % (TF.name, species.name)
        messages.info(request, msg)
        return browse_get(request)

    # use LASAGNA to align sites
    aligned, idxAligned, strands = lasagna.LASAGNA(map(lambda s: str(s.seq.lower()), site_curation_dict.keys()), 0)
    trimmed = lasagna.TrimAlignment(aligned)
    trimmed = [s.upper() for s in trimmed]
    # create weblogo for the list of sites
    weblogo_data = weblogo_uri(trimmed)

    result_dict = {'site_curation_dict':site_curation_dict,
                   'site_regulation_dict':site_regulation_dict,
                   'TF':TF,
                   'sp':species,
                   'weblogo_image_data': weblogo_data,
                   'aligned_sites': trimmed,
                   }
    response_dict = dict(make_browse_response_dict().items() + result_dict.items())
                        
    return render(request, "browse_results.html", response_dict,
                  context_instance=RequestContext(request))
    

def weblogo_uri(sequences):
    """Generate the weblogo and make it ready for direct embed into response html"""
    image_data = weblogo(sequences)
    encoded = b64encode(image_data)
    mime = "image/png"
    return ("data:" + mime + ';' + "base64," + encoded)

