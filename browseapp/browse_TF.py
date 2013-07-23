from browse_base import *
# Similarly, browse hiearchy (TF-wise) is as follows
# browse_by_TF_main
#   browse_by_TF_family
#     browse_by_TF

def get_all_TFs():
    return models.TF.objects.order_by('name')

def get_all_TF_families():
    return models.TFFamily.objects.order_by('name')

def browse_by_TF_main(request):
    """Handler for browse by TF request"""
    TF_families = get_all_TF_families()
    return render(request,
                  "browse_tf_main.html",
                  {'TF_families': TF_families},
                  context_instance=RequestContext(request))

def browse_by_TF_family(request, TF_family_id):
    """Handler for browse TF family"""
    TF_family = models.TFFamily.objects.get(pk=TF_family_id)
    TFs = models.TF.objects.filter(family=TF_family)
    return render(request,
                  "browse_tf_family.html",
                  {
                      'TF_family': TF_family,
                      'TFs': TFs
                  },
                  context_instance=RequestContext(request))

def browse_by_TF(request, TF_id):
    """Handler for browse TF"""
    TF = models.TF.objects.get(TF_id=TF_id)
    # fetch species that have curation data on this TF
    csi = models.Curation_SiteInstance.objects.filter(curation__TF=TF)
    species_ids = csi.values_list('site_instance__genome__strain', flat=True).distinct()
    species = models.Strain.objects.filter(taxonomy_id__in=species_ids).order_by('name')

    num_site_instances = {} # dictionary of num_site_instances by species_id
    num_curations = {}      # dictionary of num_curations by species_id
    # get number of sites and curations for each species
    for sp in species:
        # get filtered curation_site_instances
        filtered_csi = csi.filter(site_instance__genome__strain=sp)
        num_site_instances[sp.taxonomy_id] = filtered_csi.values_list('site_instance', flat=True).distinct().count()
        num_curations[sp.taxonomy_id] = filtered_csi.values_list('curation', flat=True).distinct().count()

    return render(request,
                  "browse_tf.html",
                  {
                      'TF': TF,
                      'species': species,
                      'num_site_instances': num_site_instances,
                      'num_curations': num_curations
                  },
                  context_instance=RequestContext(request))

