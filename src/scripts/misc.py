from core import models
from core import entrez_utils


def curations_by_date():
    """Prints all curation IDs and creation dates."""
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


def batch_assign_GO_terms_to_TF_instances(TF_name, GO_term_id):
    """Assigns given GO term to all TF instances of the given TF."""
    TF = models.TF.objects.get(name=TF_name)
    GO_term = models.GeneOntologyTerm.objects.get(GO_term_id=GO_term_id)

    for TF_instance in models.TFInstance.objects.filter(TF=TF):
        print TF_instance
        TF_instance.GO_term = GO_term
        TF_instance.save()


def check_curations_with_external_db():
    """Checks which curations have links to external databases."""
    curations = models.Curation.objects.exclude(curation_externaldatabase__isnull=True)
    print curations


def clear_spaces_on_pmids():
    """Removes trailing whitespaces on PMIDs."""
    for pub in models.Publication.objects.all():
        if pub.pmid:
            pub.pmid = pub.pmid.replace(' ', '')
            pub.save()


def batch_validate_curations():
    """Validates any curations by superusers"""
    curations = models.Curation.objects.filter(
        curator__user__username__in=['ivanerill', 'dinara'])
    sefa = models.Curator.objects.get(user__username='sefa1')
    for curation in curations:
        if not curation.validated_by:
            print curation
        curation.NCBI_submission_ready = True
        curation.validated_by = sefa
        curation.save()



def run():
    #batch_assign_GO_terms_to_TF_instances('LexA', 'GO:0006974')
    #batch_assign_GO_terms_to_TF_instances('Fur', 'GO:0071281')
    #check_curations_with_external_db()
    #clear_spaces_on_pmids()
    batch_validate_curations()
