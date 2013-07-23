from browse_base import *

# Browse hierarchy (species-wise) is as follows
# browse_by_species_main
#   browse_by_taxon
#     browse_by_species

def get_all_strains():
    return models.Strain.objects.order_by('name')

def browse_by_species_main(request):
    """Handler for browse by species. Return the taxonomy."""
    strains = get_all_strains()
    taxon = bioutils.get_all_taxon_info([strain.taxonomy_id for strain in strains])
    return render(request,
                  "browse_species_main.html",
                  {'taxon_names': sorted(list(set(taxon.values())))},
                  context_instance=RequestContext(request))

def browse_by_species_taxon(request, taxon_name):
    """Handler for browse by taxonomy. Return species only with having <taxon_name> in
    their higher level hierarchy"""
    strains = get_all_strains()
    taxon = bioutils.get_all_taxon_info([strain.taxonomy_id for strain in strains])
    filtered_strains = [strain for strain in strains if taxon[strain.taxonomy_id]==taxon_name]
    return render(request,
                  "browse_species_taxon.html",
                  {'taxon_elm': taxon_name,
                   'filtered_strains': filtered_strains
                  },
                  context_instance=RequestContext(request))

def browse_by_species(request, sp_tax_id):
    """Handler for browse by species. For the selected species (indicated by
    sp_tax_id), return the list of TFs and link to the corresponding result page."""
    sp = models.Strain.objects.get(pk=sp_tax_id)
    # fetch TFs that have curation data with this strain
    csi = models.Curation_SiteInstance.objects.filter(curation__site_instances__genome__strain=sp)
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
                  {'taxon_name': bioutils.get_taxon_info_from_file(sp_tax_id),
                   'sp': sp,
                   'TFs': TFs,
                   'num_site_instances': num_site_instances,
                   'num_curations': num_curations,
                  },
                  context_instance=RequestContext(request))

