from browse_base import *
import Queue
import sys
import time

def get_TFs(TF_families):
    """given list of families, return TFs that belong to these families"""
    return models.TF.objects.filter(family__in=TF_families)

def get_species(taxons):
    """ Given list of taxonomy elments, return all species that are children of any
    elm in taxons
    """
    species = []
    # get all leaves
    Q = Queue.Queue()
    [Q.put(t) for t in taxons]
    while not Q.empty():
        t = Q.get()
        children = t.taxonomy_set.all()
        if children:
            [Q.put(child) for child in children]
        else: # leave
            species.append(t)
    return species

def browse_TF_and_species_selected_non_motif(request, TF_param, TF_ids, species_param, species_ids):
    return browse_TF_and_species_selected(request, TF_param, TF_ids, species_param, species_ids, integrate_non_motif=True)

def browse_TF_and_species_selected(request, TF_param, TF_ids, species_param, species_ids, integrate_non_motif=False):
    """Given a list of TF/TF families and species/other taxonomy levels, return query result"""
    # get TFs
    TF_ids = TF_ids.split(',')
    if TF_param == 'TF':
        TFs = models.TF.objects.filter(TF_id__in=TF_ids)
    elif TF_param == 'TF_family':
        TF_families = models.TFFamily.objects.filter(TF_family_id__in=TF_ids)
        TFs = get_TFs(TF_families)
    # get Species
    species_ids = species_ids.split(',')
    if species_param=='species':
        species = models.Taxonomy.objects.filter(pk__in=species_ids)
    elif species_param=='taxons':
        taxons = models.Taxonomy.objects.filter(pk__in=species_ids)
        species = get_species(taxons)
            
    assert species and TFs
    # at this point we have
    # - TFs (list of TF objects)
    # - species (list of Taxonomy objects of rank='species')

    # better to get all in one query (assuming most combinations return empty)
    csis = models.Curation_SiteInstance.objects.filter(site_instance__genome__taxonomy__in=species,
                                                       curation__TF__in=TFs,
                                                       is_motif_associated=True)
    # get all non motif curation-site instances
    non_motif_csis = None
    if integrate_non_motif:
        non_motif_csis = models.Curation_SiteInstance.objects.filter(site_instance__genome__taxonomy__in=species,
                                                                     curation__TF__in=TFs,
                                                                     is_motif_associated=False)

    values = csis.values('curation__TF',
                         'site_instance__genome__taxonomy').distinct()
    reports = []
    for val in values:
        print val
        TF = val['curation__TF']
        species = val['site_instance__genome__taxonomy']
        
        filtered_csis = csis.filter(site_instance__genome__taxonomy=species,
                                    curation__TF=TF)
        filtered_non_motif_csis = None
        if non_motif_csis:
            filtered_non_motif_csis = non_motif_csis.filter(site_instance__genome__taxonomy=species,
                                                            curation__TF=TF)
        if filtered_csis:
            report = get_sites_by_TF_and_species(filtered_csis, filtered_non_motif_csis)
            report['TF'] = models.TF.objects.get(pk=TF)
            report['species'] = models.Taxonomy.objects.get(pk=species)
            reports.append(report)
            
    reports.sort(key=lambda x: x['TF'].name)
    
    #create ensemble report
    ensemble_meta_sites = []
    for report in reports:
        ensemble_meta_sites.extend(report['meta_sites'].values())
    # lasagna alignment for ensemble
    aligned, idxAligned, strands = lasagna.LASAGNA(map(lambda s:str(s[0].site_instance.seq).lower(), ensemble_meta_sites), 0)
    trimmed = lasagna.TrimAlignment(aligned) if len(aligned) > 1 else aligned
    trimmed = [s.upper() for s in trimmed]
    # create weblogo for the list of sites
    weblogo_data = bioutils.weblogo_uri(trimmed)
    
    ensemble_report = {'meta_sites': ensemble_meta_sites,
                       'aligned_sites': trimmed,
                       'weblogo_image_data': weblogo_data}

    return render(request, 'view_report.html',
                  {
                      'TF_param': TF_param,                     #
                      'TF_ids': ','.join(TF_ids),               # # # Required for redirect from/to non-motif integration page
                      'species_param': species_param,           # # 
                      'species_ids': ','.join(species_ids),     #
                      'integrate_non_motif': integrate_non_motif,
                      'reports': reports,
                      'ensemble_report': ensemble_report
                  },
                  context_instance=RequestContext(request))
    


