"""This file contains the view function for browsing by taxonomy"""

from django.shortcuts import render_to_response
from django.template import RequestContext
import browse.models as models
import browse.motif_report as motif_report
from django.shortcuts import get_object_or_404

def browse_taxonomy(request):
    """View function for browse by taxonomy. Returns the taxonomy"""
    # get all phyla
    taxonomy = dict(phyla=models.Taxonomy.objects.filter(rank='phylum')\
                    .order_by('name'))
    return render_to_response('browse_tax.html', {'taxonomy': taxonomy},
                              context_instance=RequestContext(request))

def get_species(tax_type, taxid):
    """Given a taxonomy id with its type (all, group or species), return the
    list of species"""
    if tax_type == 'all':
        return models.Taxonomy.objects.all()
    elif tax_type == 'group' or tax_type == 'species':
        tax = get_object_or_404(models.Taxonomy, pk=taxid)
        return tax.get_all_species()

def get_results_taxonomy(request, taxid):
    """Given a taxonomy id, find all species under that taxon, and return
    results for those species"""

    tax = get_object_or_404(models.Taxonomy, pk=taxid)
    # get all species
    all_species = tax.get_all_species()
    # get all curation-site-instance objects for browsed taxon
    cur_site_insts = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in=all_species)
    # generate all reports
    reports = motif_report.make_reports(cur_site_insts)

    return render_to_response("browse_results.html",
                              dict(title=tax.name,
                                   description='',
                                   reports=[report.generate_browse_result_dict()
                                            for report in reports],
                                   tax_param_type='group',
                                   tax_param=taxid,
                                   tf_param_type='all',
                                   tf_param=-1,
                                   tech_param_type='all',
                                   tech_param=-1),
                              context_instance=RequestContext(request))

