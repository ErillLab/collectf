from core import models
from core import entrez_utils

def curations_by_date():
    for curation in models.Curation.objects.all():
        date = str(curation.created).split()[0]
        print '%d,%s' % (curation.pk, date)


def fill_gene_type_field():
    """Populates the gene_type field for all genes."""
    for genome in models.Genome.objects.all():
        print genome.genome_accession
        rec = entrez_utils.get_genome(genome.genome_accession)
        genes = entrez_utils.get_genes(rec)
        gene_dict = {gene['locus_tag']: gene['gene_type']
                     for gene in genes}
        for gene in models.Gene.objects.filter(genome=genome).all():
            if gene.locus_tag in gene_dict:
                gene.gene_type = gene_dict[gene.locus_tag]
                gene.save()
            else:
                print "Missing locus tag", gene.locus_tag


def run():
    #curations_by_date()
    fill_gene_type_field()
