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


""" Function to return rendered list_all_curations page
    Function will be called by Django after user types/links to collectf.org/browse/list_all_curations
    This will first prompt Django to disambiguate the call, looking at src/browse/urls.py (the URL mapper for the browse app),
    and point it to the browse.list_all.list_all_curations function (within the browse app, in the list_all module).
    The function then performs a SQL query on the database (models is a module of the core app), and stores into curations variable
    the result of the query.
    It then returns with a request to render the 'list_all_curations.html' template (in src/templates/browse), passing it the
    curation object containing the list of curations retrieved from the database.
"""
def list_all_curations(request):
    """View function to see all curations."""
    curations = models.Curation.objects.all().order_by('-curation_id')
    return render(request, 'list_all_curations.html', {"curations": curations})
