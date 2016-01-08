"""The view functions for browsing by TF and TF family."""

from django.shortcuts import render
from django.shortcuts import get_object_or_404

from core import models
from .browse_by_utils import curation_site_instances_values_list


def browse_TF(request):
    """Returns the TF treeview for browsing by TF and family."""
    context = {'TF_families': models.TFFamily.objects.all().order_by('name')}
    return render(request, 'browse_by_TF.html', context)


def get_results_by_TF_family(request, object_id):
    """Returns TF-species pairs that have binding sites for a given family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    TF_species_pairs = curation_site_instances_values_list(
        models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF__in=TF_family.TFs.all()))
    return render(request, 'browse_results.html',
                  {'title': TF_family.name,
                   'description': TF_family.description,
                   'TF_species_pairs': TF_species_pairs})


def get_results_by_TF(request, object_id):
    """Returns TF-species pairs that have binding sites for a given TF."""
    TF = get_object_or_404(models.TF, TF_id=object_id)
    TF_species_pairs = curation_site_instances_values_list(
        models.Curation_SiteInstance.objects.filter(
            curation__TF_instances__TF=TF))
    return render(request, 'browse_results.html',
                  {'title': TF.name,
                   'description': TF.description,
                   'TF_species_pairs': TF_species_pairs})
