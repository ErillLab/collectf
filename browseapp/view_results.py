import collectfapp.views
from browse_base import *
from django.core.urlresolvers import reverse
from django.http import HttpResponseRedirect
from collections import namedtuple

def prepare_results(motif_csi_list, non_motif_csi_list, integrate_non_motif=False):
    """
    Given lists of motif-associated and non-motif-associated curation site
    instances, pass them to the view results function.
    """
    csis = models.Curation_SiteInstance.objects.filter(pk__in=motif_csi_list)
    # get non-motif-associated curation-site-instances
    # if non_motif_csi_list is empty filter doesn't work!
    non_motif_csis = models.Curation_SiteInstance.objects.none()
    if integrate_non_motif and non_motif_csi_list:
        non_motif_csis = models.Curation_SiteInstance.objects.filter(pk__in=non_motif_csi_list)
    # get data, associated with curation-site-instances
    values = csis.values('curation__TF','site_instance__genome__taxonomy').distinct()
    # create all reports
    reports = []
    for val in values:
        TF = val['curation__TF']
        species = val['site_instance__genome__taxonomy']
        filtered_csis = csis.filter(site_instance__genome__taxonomy=species, curation__TF=TF)
        if not filtered_csis: continue
        filtered_non_motif_csis = non_motif_csis.filter(site_instance__genome__taxonomy=species, curation__TF=TF)
        meta_sites = get_sites_by_TF_and_species(filtered_csis,filtered_non_motif_csis)
        report = meta_sites_to_report(meta_sites)
        report['TF'] = models.TF.objects.get(pk=TF)
        report['species'] = models.Taxonomy.objects.get(pk=species)
        reports.append(report)
    reports.sort(key=lambda x: x['TF'].name)
    # ensemble report
    ensemble_meta_sites = [meta_site['sites'][0] for report in reports
                           for meta_site in report['meta_sites']]
    # lasagna alignment for ensemble
    trimmed = bioutils.call_lasagna(map(lambda s: s.site_instance, ensemble_meta_sites))
    # create weblogo for the list of sites
    ensemble_report = {'unaligned_sites': map(lambda s: s.site_instance.seq, ensemble_meta_sites),
                       'aligned_sites': trimmed}
    return {'reports': reports,
            'ensemble_report': ensemble_report,
            'motif_csi_list': ','.join(map(str, motif_csi_list)),
            'non_motif_csi_list': ','.join(map(str, non_motif_csi_list)) if non_motif_csi_list else ""}
    
def view_results(request):
    #assert request.POST, '
    if not request.POST:
        #no GET handler for this function
        return HttpResponseRedirect(reverse(collectfapp.views.home))
        
    csi_list = request.POST['motif_csi_list'].strip().split(',')
    non_motif_csi_list = None
    if request.POST['non_motif_csi_list']:
        non_motif_csi_list = request.POST['non_motif_csi_list'].strip().split(',')
    integrate_non_motif = bool('integrate_non_motif' in request.POST)
    template = prepare_results(csi_list, non_motif_csi_list, integrate_non_motif)
    return render(request,
                  'view_report.html',
                  template,
                  context_instance=RequestContext(request))

def get_sites_by_TF_and_species(curation_site_instances, non_motif_curation_site_instances):
    if not curation_site_instances: return None
    meta_sites = group_curation_site_instances(curation_site_instances, non_motif_curation_site_instances)
    return meta_sites

def meta_sites_to_report(meta_sites):
    """Given a list of meta-site objects, return dictionary in proper format"""
    # Use LASAGNA to align sites
    trimmed = bioutils.call_lasagna(map(lambda ms: ms['sites'][0].site_instance, meta_sites))
    return {
        'meta_sites': meta_sites,
        'aligned_sites': trimmed,
    }

def group_curation_site_instances(curation_site_instances, non_motif_curation_site_instances):
    """Given a queryset of motif-associated and non-motif-associated
    curation-site-instance objects, group them by meta-site-instances as well as the
    regulated genes"""
    # meta_site is dictionary of {sites, curations, regulations, techniques}
    def group_motif_associated_site_instances():
        meta_sites = []
        # group motif-associated curation-site-instances
        for csi in curation_site_instances.iterator():
            for ms in meta_sites: # check all existing meta-sites if there is any fit
                if (ms['sites'][0].site_instance.genome == csi.site_instance.genome and
                    ms['sites'][0].curation.TF_instance == csi.curation.TF_instance and
                    bioutils.overlap_site_meta_site(csi, ms['sites'])):
                    ms['sites'].append(csi)
                    break
            else: # that means it doesn't overlap with any of the existing metasites
                new_meta_site = dict(sites=[csi], curations=[], regulations=[], techniques=[])
                meta_sites.append(new_meta_site) # new meta-site with one curation-site-instance
        return meta_sites
    
    def integrate_non_motif_associated_site_instances(meta_sites):
        # integrate non-motif-associated curation-site-instances
        if not non_motif_curation_site_instances: return
        for ncsi in non_motif_curation_site_instances.iterator():
            for ms in meta_sites:
                if (ms['sites'][0].site_instance.genome == ncsi.site_instance.genome and
                    ms['sites'][0].curation.TF_instance == ncsi.curation.TF_instance and
                    bioutils.overlap_non_motif_site_meta_site(ncsi, ms['sites'])):
                    ms['sites'].append(ncsi)
                    # no break here; incorporate each non-motif-associated site to all
                    # possible meta-sites
            else: pass # nothing to do

    def group_curations_for_each_meta_site(meta_sites):
        for ms in meta_sites:
            ms['curations'] = list(set([csi.curation for csi in ms['sites']]))

    def group_techniques_for_each_meta_site(meta_sites):
        all_techniques = models.ExperimentalTechnique.objects.all()
        for ms in meta_sites:
            ms['techniques'] = list(set(all_techniques.filter(preset_function__in=['binding', 'expression'],
                                                              curation__in=ms['curations'])))

    def group_regulations_for_each_meta_site(meta_sites):
        #all_regs = models.Regulation.objects.all()
        all_exp_verified_regs = models.Regulation.objects.filter(evidence_type="exp_verified")
        all_inferred_regs = models.Regulation.objects.filter(evidence_type="inferred")
        for ms in meta_sites:
            regs_exp_verified = all_exp_verified_regs.filter(curation_site_instance__in=ms['sites'])
            regs_inferred = all_inferred_regs.filter(curation_site_instance__in=ms['sites'])
            genes_exp_verified = regs_exp_verified.values_list('gene', flat=True)
            genes_inferred = regs_inferred.values_list('gene', flat=True)
            # filter multiple regulations of the same gene
            ms['regulations'] = [reg for i,(gene,reg) in enumerate(zip(genes_exp_verified, regs_exp_verified))
                                 if gene not in genes_exp_verified[:i]]
            ms['regulations'].extend([reg for i,(gene,reg) in enumerate(zip(genes_inferred, regs_inferred))
                                      if gene not in genes_inferred[:i] and gene not in genes_exp_verified])

    meta_sites = group_motif_associated_site_instances()
    integrate_non_motif_associated_site_instances(meta_sites)
    group_curations_for_each_meta_site(meta_sites)
    group_techniques_for_each_meta_site(meta_sites)
    group_regulations_for_each_meta_site(meta_sites)
    return meta_sites
