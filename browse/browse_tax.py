"""
The view function for browsing by taxonomy.
"""

from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404
from django.shortcuts import render_to_response
from django.template import RequestContext

import browse.models as models
import browse.motif_report as motif_report
from browse.static_reports import get_static_reports
from browse.view_reports import view_reports_by_taxonomy

def browse_taxonomy(request):
    """View function for browse by taxonomy. Returns the taxonomy"""
    # get all phyla
    taxonomy = {'phyla': models.Taxonomy.objects.filter(rank='phylum')}
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

def get_results_taxonomy(request, object_id):
    """Returns motif reports for a given taxonomy ID."""
    tax = get_object_or_404(models.Taxonomy, taxonomy_id=object_id)
    reports, _ = get_static_reports('taxonomy_%s' % object_id)
    return render_to_response(
        'browse_results.html',
        {'title': tax.name,
         'description': '',
         'reports': [r.generate_browse_result_dict() for r in reports],
         'combined_report_url': reverse(view_reports_by_taxonomy,
                                        args=(tax.taxonomy_id,))},
        context_instance=RequestContext(request))

