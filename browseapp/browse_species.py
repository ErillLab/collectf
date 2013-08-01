from browse_base import *

# Browse hierarchy (species-wise) is as follows
# browse_by_species_main
#   browse_by_taxon
#     browse_by_species

def get_all_species():
    species_ids = models.Genome.objects.values_list('taxonomy', flat=True).distinct()
    return models.Taxonomy.objects.filter(pk__in=species_ids)

def browse_by_species_main(request):
    """Handler for browse by species. Return the taxonomy."""
    species = get_all_species()
    orders = sorted(list(set(sp.get_order() for sp in species)), key=lambda x:x.name)
    
    return render(request,
                  "browse_species_main.html",
                  {'species': orders},
                  context_instance=RequestContext(request))

def browse_by_species_taxon(request, tax_id):
    """Handler for browse by taxonomy. Return species only with having <taxon_id> in
    their higher level hierarchy"""
    taxon = models.Taxonomy.objects.get(pk=tax_id)
    species = get_all_species().filter(parent__parent__parent=tax_id)
    return render(request,
                  "browse_species_taxon.html",
                  {'taxon': taxon,
                   'species': species
                  },
                  context_instance=RequestContext(request))

def browse_by_species(request, sp_tax_id):
    """Handler for browse by species. For the selected species (indicated by
    sp_tax_id), return the list of TFs and link to the corresponding result page."""
    sp = models.Taxonomy.objects.get(pk=sp_tax_id)
    # fetch TFs that have curation data with this strain
    csi = models.Curation_SiteInstance.objects.filter(curation__site_instances__genome__taxonomy=sp)
    TF_ids = csi.values_list('curation__TF', flat=True).distinct()
    TFs = models.TF.objects.filter(TF_id__in=TF_ids).order_by('name')
    num_site_instances= {} # dictionary of num_site_instances by TF_instance_id
    num_curations = {}     # dictionary of num_curations by TF_instance_id
    # get number of sites and curations for each TF
    for TF in TFs:
        filtered_csi = csi.filter(curation__TF=TF)
        num_site_instances[TF.TF_id] = filtered_csi.values_list('site_instance', flat=True).distinct().count()
        num_curations[TF.TF_id] = filtered_csi.values_list('curation', flat=True).distinct().count()

    return render(request,
                  "browse_species.html",
                  {'sp': sp,
                   'parent': sp.get_order(),
                   'TFs': TFs,
                   'num_site_instances': num_site_instances,
                   'num_curations': num_curations,
                  },
                  context_instance=RequestContext(request))

