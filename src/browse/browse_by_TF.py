"""The view functions for browsing by TF and TF family."""

from django.shortcuts import render
from django.shortcuts import get_object_or_404

from core import models

def browse_by_TF(request):
    """Returns the TF treeview for browsing by TF and family."""
    context = {'TF_families': models.TFFamily.objects.all().order_by('name')}
    print context
    return render(request, 'browse_by_TF.html', context)

def get_results_by_TF_family(request, object_id):
    """Returns TF-species pairs that have binding sites for a given family."""
    TF_family = get_object_or_404(models.TFFamily, TF_family_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF__in=TF_family.TFs.all())
    TF_species_pairs = curation_site_instances.values_list(
        'curation__TF_instances__TF__TF_id',
        'curation__TF_instances__TF__name',
        'site_instance__genome__taxonomy__taxonomy_id',
        'site_instance__genome__taxonomy__name').distinct()
    return render(request, 'browse_results.html',
                  {'title': TF_family.name,
                   'description': TF_family.description,
                   'TF_species_pairs': TF_species_pairs})

def get_results_by_TF(request, object_id):
    """Returns TF-species pairs that have binding sites for a given TF."""
    TF = get_object_or_404(models.TF, TF_id=object_id)
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        curation__TF_instances__TF=TF)
    TF_species_pairs = curation_site_instances.values_list(
        'curation__TF_instances__TF__TF_id',
        'curation__TF_instances__TF__name',
        'site_instance__genome__taxonomy__taxonomy_id',
        'site_instance__genome__taxonomy__name').distinct()
    return render(request, 'browse_results.html',
                  {'title': TF.name,
                   'description': TF.description,
                   'TF_species_pairs': TF_species_pairs})

