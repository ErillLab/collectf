"""This file contains the view function for browsing by taxonomy"""

from django.shortcuts import render
from django.shortcuts import render_to_response
from django.http import HttpResponse
from django.http import HttpResponseRedirect
from django.template import RequestContext
import models
import Queue
import motif_report

def browse_taxonomy(request):
    """View function for browse by taxonomy. Returns the taxonomy"""
    # get all phyla
    taxonomy = {'phyla': models.Taxonomy.objects.filter(rank='phylum').order_by('name')}
    return render_to_response('browse_tax.html', {'taxonomy': taxonomy},
                              context_instance=RequestContext(request))

def get_results_taxonomy(request, taxid):
    """Given a taxonomy id, find all species under that taxon, and return
    results for those species"""

    def get_all_sp(taxid):
        """Given a taxon t, return all species (leaves) of the phylogeny where t
        is the root of the tree."""
        all_sp = []
        Q = Queue.Queue()
        Q.put(models.Taxonomy.objects.get(pk=taxid))
        while not Q.empty():
            node = Q.get()
            children =  node.taxonomy_set.all()
            if children:
                for c in children:
                    Q.put(c)
            else:
                all_sp.append(node)
        return all_sp

    # get all species
    all_species = get_all_sp(taxid)
    # get all curation-site-instance objects for browsed taxon
    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in = all_species)

    # generate all reports
    reports = motif_report.make_reports(cur_site_insts)

    return render_to_response("browse_results.html",
                              {'title': models.Taxonomy.objects.get(pk=taxid).name,
                               'description': '',
                               'all_cur_site_insts': [pk for report in reports
                                                      for pk in report.get_all_cur_site_insts_ids()],
                               'reports': [report.generate_browse_result_dict() for report in reports],},
                                context_instance=RequestContext(request))
    
