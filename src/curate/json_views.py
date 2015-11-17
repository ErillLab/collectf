import json
from django.http import JsonResponse

from core import models


def get_genomes(request):
    """Gets the list of all genomes."""
    genomes =  list(models.Genome.objects.values(
        'genome_accession', 'organism').distinct())
    return JsonResponse(genomes, safe=False)


def get_TF_instances(request):
    """Gets the list of all TF instances."""
    TF_instances = list(models.TFInstance.objects.values(
        'uniprot_accession', 'refseq_accession', 'description').distinct())
    return JsonResponse(TF_instances, safe=False)


def uniprot_to_refseq(request, uniprot_accession):
    """Returns the RefSeq accession for the given UniProt accession"""
    # Check if TF-instance is in database.
    refseq = ''
    try:
        TF_instance = models.TFInstance.objects.get(
            uniprot_accession=uniprot_accession)
        refseq = TF_instance.refseq_accession
    except models.TFInstance.DoesNotExist:
        pass

    # TODO(sefa): Try to fetch using UniProt mapping service.
    
    return JsonResponse(refseq, safe=False)
