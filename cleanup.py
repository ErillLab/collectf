
import sys
from collectfapp.models import *
import time


def get_all_genes():
    return Gene.objects.all()

def get_all_regulations():
    return Regulation.objects.all()

def regulation_duplicates():

    all_genes = get_all_genes()
    all_regulations = get_all_regulations()
    for reg in regulations:
        same_genes = [g for g in all_genes if g.gene_accession==reg.gene.gene_accession]
        print same_genes
