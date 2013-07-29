from browse_base import *

def browse_TF_and_species_selected(request, TF_id, species_id):
    TF = models.TF.objects.get(TF_id=TF_id)
    species = models.Strain.objects.get(pk=species_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__strain=species,
        curation__TF=TF,
        is_motif_associated=True) # get motif_associated ones for now.
    return get_sites_by_TF_and_species(request,TF,species,curation_site_instances)

def browse_TF_and_species(request):
    return browse_TF_and_species_get(request) if not request.POST else browse_TF_and_species_post(request)

def browse_TF_and_species_get(request):
    """Handler for TF/species browse page. Renders the page with the form which
    enables user to query the database."""
    return render(request, "browse.html", get_template_dict(), context_instance=RequestContext(request))

def browse_TF_and_species_post(request):
    """Process form to get selected TF, species and experimental techniques. Retrieve
    binding sites for selected TF and species from database, filter them by
    experimental techniques, return to the user."""
    
    TF = models.TF.objects.get(TF_id=request.POST['TF'])
    species = models.Strain.objects.get(pk=request.POST['species'])
    experimental_techniques_1 = request.POST['tech1']
    experimental_techniques_2 = request.POST['tech2']
    experimental_techniques_3 = request.POST['tech3']
    boolean1 = request.POST['boolean1']
    boolean2 = request.POST['boolean2']
    assert boolean1 in ['and', 'or'] and boolean2 in ['and', 'or']
    assert 'tech1' in request.POST and 'tech2' in request.POST and 'tech3' in request.POST

    q1 = techniques_JSON_to_Q(experimental_techniques_1)
    q2 = techniques_JSON_to_Q(experimental_techniques_2)
    q3 = techniques_JSON_to_Q(experimental_techniques_3)

    #  get curation objects
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF=TF,
        site_instance__genome__strain=species,
        is_motif_associated=True)
    # filter them by experimental techniques
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
    return get_sites_by_TF_and_species(request, TF, species, curation_site_instances)

def get_sites_by_TF_and_species(request, TF, species, curation_site_instances):
    meta_sites, meta_site_curation_dict, meta_site_regulation_dict = group_curation_site_instances(curation_site_instances)
    # if there is no site, message
    if not meta_sites:
        msg = "No site found for transcription factor %s in the genome of %s." % (TF.name, species.name)
        messages.info(request, msg)
        return browse_TF_and_species_get(request)

    # Use LASAGNA to align sites    # use LASAGNA to align sites
    aligned, idxAligned, strands = lasagna.LASAGNA(map(lambda s:str(s[0].site_instance.seq).lower(), meta_sites.values()), 0)
    trimmed = lasagna.TrimAlignment(aligned) if len(aligned) > 1 else aligned
    print trimmed
    trimmed = [s.upper() for s in trimmed]
    # create weblogo for the list of sites
    weblogo_data = bioutils.weblogo_uri(trimmed)

    result_dict = get_template_dict()
    result_dict['meta_sites'] = meta_sites
    result_dict['meta_site_curation_dict'] = meta_site_curation_dict
    result_dict['meta_site_regulation_dict'] = meta_site_regulation_dict
    result_dict['TF'] = TF
    result_dict['species'] = species
    result_dict['weblogo_image_data'] = weblogo_data
    result_dict['aligned_sites']= trimmed
    return render(request, "browse_results.html", result_dict, context_instance=RequestContext(request))

def overlap_site_meta_site(site_instance, meta_site_instance, overlap_th=0.8):
    """Given a site instance (site_instance) and a list of site instances
    (meta_site_instance), return whether site instance overlaps enough with any site
    instance in the meta_site_instance."""
    return any(bioutils.get_overlap((site_instance.start, site_instance.end),
                                    (msi.start, msi.end)) for msi in meta_site_instance)

    
def group_curation_site_instances(curation_site_instances):
    # Group curation_site_instance objects by meta site-instances
    # group all curation_site_instances
    meta_sites = dict()
    for csi in curation_site_instances:
        # search for a meta-site-instance
        for i in meta_sites.keys():
            if overlap_site_meta_site(csi.site_instance, [m.site_instance for m in meta_sites[i]]):
                meta_sites[i].append(csi)
                break
        else:
            meta_sites[len(meta_sites)+1] = [csi]
    # 
    meta_site_curation_dict = {}
    meta_site_regulation_dict = {}

    for ms_id in meta_sites:
        meta_site_curation_dict[ms_id] = []
        for csi in meta_sites[ms_id]:
            meta_site_curation_dict[ms_id].append(csi.curation)
            meta_site_regulation_dict[ms_id] = []
            for reg in csi.regulation_set.all():
                # check if regulated gene is already in the list (from another curation)
                same_reg_gene = [r for r in meta_site_regulation_dict[ms_id] if r.gene==reg.gene]
                if same_reg_gene and same_reg_gene.evidence_type=='inferred':
                    same_reg_gene[0].evidence_type = reg.evidence_type
                if not same_reg_gene:
                    meta_site_regulation_dict[ms_id].append(reg)
                    
    return meta_sites, meta_site_curation_dict, meta_site_regulation_dict

def get_template_dict():
    """Pass dictionary to the browse HTML template page."""
    # group techniques
    binding_techniques = {}
    expression_techniques = {}
    insilico_techniques = {}

    all_categories = models.ExperimentalTechniqueCategory.objects.all()
    for category in all_categories:
        # find all techniques that belong to that category
        # category and techniques have n:n relationship
        techs = models.ExperimentalTechnique.objects.filter(categories=category)
        binding_techniques[category.name] = techs.filter(preset_function='binding')
        expression_techniques[category.name] = techs.filter(preset_function='expression')
        insilico_techniques[category.name] = techs.filter(preset_function='insilico')

    # remove empty keys from dict
    binding_techniques = dict((x,y) for (x,y) in binding_techniques.items() if y)
    expression_techniques = dict((x,y) for (x,y) in expression_techniques.items() if y)
    insilico_techniques = dict((x,y) for (x,y) in insilico_techniques.items() if y)
    
    return dict(TFs=models.TF.objects.all().order_by('name'),
                species=models.Strain.objects.all().order_by('name'),
                binding_techniques=binding_techniques,
                expression_techniques=expression_techniques,
                insilico_techniques=insilico_techniques)

def techniques_JSON_to_Q(JSON_string):
    """Given JSON string received from the form, parse the techniques (names), get
    the ids for them, build the Q object for filtering"""
    j = json.loads(JSON_string)
    techniques = map(lambda x: map(lambda y: y['key'], x['values']), j)
    techniques = [t for grp in techniques for t in grp] # flatten list
    q = Q(curation__curation_id=-9999)
    for t in techniques:
        technique_object = models.ExperimentalTechnique.objects.get(name=t)
        q = q | Q(curation__experimental_techniques=technique_object)
    return q

def export_sites(request):
    """Given a list of sites, report FASTA/CSV file
    containing sites for particular TF and species"""
    if 'fasta' in request.POST: export_format = 'fasta'
    elif 'csv'in request.POST: export_format = 'csv'
    assert export_format
    site_ids = request.POST.getlist('site_id')
    sites = models.MetaSiteInstance.objects.filter(pk__in=site_ids)
    filename = "sites.fasta" if export_format == 'fasta' else "sites.csv"
    # set HttpResponse stuff
    response = HttpResponse(content_type='application/download')
    response['Content-Disposition'] = 'attachment;filename="%s"' % filename
    # write all sites to file
    if export_format == 'csv':
        response.write('\t'.join(['genome', 'site_start', 'site_end', 'site_seq', 'regulated_genes', 'references']))
        response.write('\n')
    for site in sites:
        response.write(site.to_fasta() if export_format=='fasta' else site.to_csv())
        if export_format == 'csv':
            regulations = models.Regulation.objects.filter(curation_site_instance__meta_site_instance=site)
            response.write('\t' + ','.join(reg.gene.name for reg in regulations))
            response.write('\n')
    return response
