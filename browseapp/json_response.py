from browse_base import *

def get_genomes(request):
    genomes = list(models.Genome.objects.values('genome_accession', 'organism').distinct())
    return HttpResponse(json.dumps(genomes), mimetype="application/json")

def get_TF_instances(request):
    TF_instances = list(models.TFInstance.objects.values('protein_accession', 'name', 'description').distinct())
    return HttpResponse(json.dumps(TF_instances), mimetype="application/json")
