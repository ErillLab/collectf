import models
import json
from django.http import HttpResponse

def get_genomes(request):
    """Get the list of all genomes"""
    genomes = list(
        models.Genome.objects.values('genome_accession', 'organism').distinct())
    return HttpResponse(json.dumps(genomes), mimetype="application/json")

def get_TF_instances(request):
    """Get the list of all TF instances"""
    TF_instances = list(models.TFInstance.objects.values(
            'uniprot_accession', 'refseq_accession', 'description').distinct())
    return HttpResponse(json.dumps(TF_instances), mimetype="application/json")
