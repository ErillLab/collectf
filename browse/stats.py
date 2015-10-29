import os
from django.shortcuts import render
from django.contrib import messages
from django.core.urlresolvers import reverse
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from collectf import settings
import models
import pickle
from browse.motif_report import make_reports

DBSTATS_PICKLE_FILE = os.path.join(settings.PICKLE_ROOT,
                                   "collectf_dbstats.pickle")

def curator_roster(request):
    """Get the list of curators."""
    return render(request,
                  "curator_roster.html",
                  context_instance=RequestContext(request))

def list_tfs(request):
    """Returns all TF families."""
    tf_families = models.TFFamily.objects.all()
    return render(request,
                  'all_tfs.html',
                  {'tf_families': tf_families},
                  context_instance=RequestContext(request))

def list_species(request):
    """Returns all phyla objects, therefore all species."""
    phyla = models.Taxonomy.objects.filter(rank='phylum')
    return render(request,
                  'all_species.html',
                  {'phyla': phyla},
                  context_instance=RequestContext(request))

def list_experimental_techniques(request):
    """Get the experimental techniques."""
    exp_techniques = models.ExperimentalTechnique.objects.all()
    return render(request,
                  "all_experimental_techniques.html",
                  {'techs': exp_techniques},
                  context_instance=RequestContext(request))

def release_history(request):
    """Get release history"""
    return render(request, "release_history.html")

def publication_complete_ratio():
    """Return the percentage of publication completed."""
    return (models.Publication.objects.filter(curation_complete=True).count() *
            100.0 / models.Publication.objects.count())

def num_curations_and_sites():
    """Return the number of curations and sites for each pair of TF and
    species"""
    all_TFs = models.TF.objects.all()
    all_species = models.Taxonomy.objects.filter(rank='species')

    num_curations = {}
    num_sites = {}
    for TF in all_TFs:
        num_curations[TF.name] = {}
        num_sites[TF.name] = {}
        for species in all_species:
            cur_site_insts = models.Curation_SiteInstance.objects.filter(
                site_type__in=['motif_associated', 'var_motif_associated'],
                site_instance__genome__taxonomy=species,
                curation__TF_instances__TF=TF)
            num_cur = cur_site_insts.values('curation').distinct().count()
            num_sit = cur_site_insts.values('site_instance').distinct().count()
            num_curations[TF.name][species.name] = num_cur
            num_sites[TF.name][species.name] = num_sit
    return num_curations, num_sites

def num_sites_by_TF_species():
    """Return the number of sites for each pair of TF and species"""
    return None

def update_stats(request):
    """Update CollecTF statistics page."""
    all_TFs = models.TF.objects.all()
    all_species_ids = models.Genome.objects.values_list("taxonomy", flat=True).distinct()
    all_species = models.Taxonomy.objects.filter(pk__in=all_species_ids)
    num_curations, num_sites = num_curations_and_sites()
    d = dict(num_TFs=all_TFs.count(),
             num_species=all_species.count(),
             num_curations=models.Curation.objects.count(),
             num_sites=models.SiteInstance.objects.count(),
             num_publications=models.Publication.objects.count(),
             pub_completed='%.1f' % publication_complete_ratio(),
             num_curations_by_TF_species=num_curations,
             num_sites_by_TF_species=num_sites,
             TFs=sorted([tf.name for tf in all_TFs]),
             species=sorted([sp.name for sp in all_species]))

    pickle.dump(d, open(DBSTATS_PICKLE_FILE, 'w'))

    message = "Database statistics updated successfully."
    messages.add_message(request, messages.SUCCESS, message)
    return HttpResponseRedirect(reverse(curation_stats))

def curation_stats(request):
    """Handler for curation statistics page. Count the number of curations/sites
    for each TF and species in the database."""
    response_dict = pickle.load(open(DBSTATS_PICKLE_FILE, 'r'))
    return render(request, "database_stats.html", response_dict,
                  context_instance=RequestContext(request))

def view_all_curations(request):
    """Handler function to see all curations at once.  This function renders the
    page with the list of all curations in the database"""
    all_curations = models.Curation.objects.all().order_by('-curation_id')
    return render_to_response("view_all_curation.html",
                              {"curations": all_curations},
                              context_instance=RequestContext(request))


def view_all_publications(request):
    """Handler function to see all publications in the database.  This is for
    internal use, to see all publications in the database."""
    all_pubs = models.Publication.objects.all().order_by('-pmid')
    return render_to_response("view_all_publication.html",
                              {"publications": all_pubs},
                              context_instance=RequestContext(request))

def list_all_motifs(request):
    """Handler function to list all motifs in the database."""
    all_csis = models.Curation_SiteInstance.objects.all()
    reports = make_reports(all_csis)
    return render_to_response("list_motifs.html",
                              {'reports': reports},
                              context_instance=RequestContext(request))
