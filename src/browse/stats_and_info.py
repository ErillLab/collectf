"""View functions about CollecTF stats and info."""

from django.shortcuts import render

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


def stats(request):
    """Returns the stats page."""
    
    return render(
        request,
        'database_stats.html',
        {'TF_count': models.TF.objects.count(),
         'species_count': models.Taxonomy.objects.filter(
             rank='species').count(),
         'curation_count': models.Curation.objects.count(),
         'sites_count': models.SiteInstance.objects.count(),
         'publication_count': models.Publication.objects.count(),
         'publication_complete_ratio': publication_complete_ratio()})
