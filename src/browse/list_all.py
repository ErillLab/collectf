"""Functions to list all

- TFs
- species
- experimental techniques
- curations
- publications
"""

from django.shortcuts import render

from core import models


def list_all_TFs(request):
    """Returns the list of all TFs."""
    TF_families = models.TFFamily.objects.all()
    return render(request, 'list_all_TFs.html', {'TF_families': TF_families})


def list_all_species(request):
    """Returns the list of all species."""
    phyla = models.Taxonomy.objects.filter(rank='phylum')
    return render(request, 'list_all_species.html', {'phyla': phyla})


def list_all_experimental_techniques(request):
    """Gets the experimental techniques."""
    exp_techniques = models.ExperimentalTechnique.objects.all()
    return render(request, 'list_all_experimental_techniques.html',
                  {'techniques': exp_techniques})


def list_all_publications(request):
    """View function to see all publications in the database."""

    publications = models.Publication.objects.all().order_by('-pmid')
    return render(request, 'list_all_publications.html',
                  {"publications": publications})


def list_all_curations(request):
    """View function to see all curations."""
    curations = models.Curation.objects.all().order_by('-curation_id')
    return render(request, 'list_all_curations.html', {"curations": curations})
