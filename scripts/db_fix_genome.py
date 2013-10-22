"""Script to move genome sequences into a separate table in the database"""
from collectfapp import models
def db_fix_genome():
    for g in models.Genome.objects.all():
        print g.genome_accession
        genome_seq = models.GenomeSequence.objects.create(sequence=g.sequence)
        genome_seq.save()
        g.genome_sequence = genome_seq
        g.save()

def check():
    for g in models.Genome.objects.iterator():
        print g.genome_accession
        assert g.genome_sequence.sequence == g.sequence
        

def run():
    #db_fix_genome()
    check()
    
