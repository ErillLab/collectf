"""
The view function for browsing by taxonomy.
"""

from django.core.urlresolvers import reverse
from django.shortcuts import get_object_or_404
from django.shortcuts import render_to_response
from django.template import RequestContext

import models
from .motif_report import make_reports
from .view_reports import view_reports_by_taxonomy

def browse_taxonomy(request):
    """View function for browse by taxonomy. Returns the taxonomy"""
    # get all phyla
    taxonomy = {'phyla': models.Taxonomy.objects.filter(rank='phylum')}
    return render_to_response('browse_tax.html', {'taxonomy': taxonomy},
                              context_instance=RequestContext(request))

def get_results_taxonomy(request, object_id):
    """Returns motif reports for a given taxonomy ID."""
    tax = get_object_or_404(models.Taxonomy, taxonomy_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__taxonomy__in=tax.get_all_species())
    reports = make_reports(curation_site_instances)
    return render_to_response(
        'browse_results.html',
        {'title': tax.name,
         'description': '',
         'reports': reports,
         'combined_report_url': reverse(view_reports_by_taxonomy,
                                        args=(tax.taxonomy_id,))},
        context_instance=RequestContext(request))

