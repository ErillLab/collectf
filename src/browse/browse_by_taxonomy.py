"""The view functions for browsing by taxonomy."""

from django.core.urlresolvers import reverse
from django.shortcuts import render
from django.shortcuts import get_object_or_404

from core import models
from .browse_by_utils import curation_site_instances_values_list


def browse_taxonomy(request):
    """View function for browse by taxonomy. Returns the taxonomy."""
    taxonomy = {'phyla': models.Taxonomy.objects.filter(rank='phylum')}
    return render(request, 'browse_by_taxonomy.html', {'taxonomy': taxonomy})


def get_results_by_taxonomy(request, object_id):
    """Returns motif reports for a given taxonomy ID."""
    taxon = get_object_or_404(models.Taxonomy, pk=object_id)
    TF_species_pairs = curation_site_instances_values_list(
        models.Curation_SiteInstance.objects.filter(
            site_instance__genome__taxonomy__in=taxon.get_all_species()))
    return render(request, 'browse_results.html',
                  {'title': taxon.name,
                   'description': '',  # Will be populated from Wikipedia.
                   'TF_species_pairs': TF_species_pairs,
                   'ensemble_report_url': reverse(
                       'view_motif_reports_by_taxonomy', args=(taxon.pk,))})
