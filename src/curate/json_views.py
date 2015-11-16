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
