"""View functions about CollecTF stats and info."""
from collections import defaultdict

from django.shortcuts import render
from django.db.models import Count

from core import models


def curator_roster(request):
    """Returns the list of curators."""
    return render(request, 'curator_roster.html')


def release_history(request):
    """Returns the release history."""
    return render(request, 'release_history.html')


def publication_complete_ratio():
    """Returns the percentage of publication completed."""
    num_papers = models.Publication.objects.count()
    num_completed_papers = models.Publication.objects.filter(
        curation_complete=True).count()
    return int(100.0 * num_completed_papers / num_papers)


def curation_counts():
    """Returns the number of curations by TF and species."""
    return models.Curation_SiteInstance.objects.all().values_list(
        'site_instance__genome__organism',
        'curation__TF_instances__TF__name').annotate(n=Count('curation'))

def site_instance_stats():
    """Returns the number of site instances by TF and species."""
    return models.Curation_SiteInstance.objects.all().values(
        'site_instance__genome__organism',
        'curation__TF_instances__TF').annotate(n=Count('site_instance'))


def stats(request):
    """Returns the stats page."""
    all_TFs = models.TF.objects.values_list('name', flat=True)
    all_species = models.Genome.objects.values_list('organism', flat=True)

    return render(
        request,
        'database_stats.html',
        {'curation_count': models.Curation.objects.count(),
         'sites_count': models.SiteInstance.objects.count(),
         'publication_count': models.Publication.objects.count(),
         'publication_complete_ratio': publication_complete_ratio(),
         'TFs': all_TFs,
         'species': all_species})
            
