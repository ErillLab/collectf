# List of functions to create objects in db
from models import *
import bioutils

def citation(pubrec):
    """Create citation string from publication record"""
    title = pubrec["Title"]
    authors = pubrec["AuthorList"]
    journal = pubrec["FullJournalName"]
    return '|'.join([title, ','.join(authors), journal])

def make_pub(pubrec, cd):
    """Given publication record, retrieved from NCBI db (if pubmed publication),
    and entered user data, create models.Publication object and return it.
    """
    pmid = cd.get("pmid", None)  # None if not pubmed publication
    url ="http://www.ncbi.nlm.nih.gov/pubmed?term=%s" % pmid if pmid else cd['URL']
    publication_type = "pubmed" if pmid else "nonpubmed"
    p = Publication(publication_type = publication_type,
                    pmid=pmid,
                    citation=citation(pubrec),
                    url=url,
                    pdf=None,
                    contains_promoter_data=cd["contains_promoter_data"],
                    contains_expression_data=cd["contains_expression_data"],
                    submission_notes=cd["submission_notes"],
                    curation_complete=False)
    return p

def make_genome(genome_record):
    """Given genome record that is fetched from NCBI, create a models.Genome object"""
    gc = bioutils.GC(genome_record.seq)  # GC content
    # get or create strain if not in database
    org = bioutils.get_org_name(genome_record)           # organism name
    strain_taxon = bioutils.get_org_taxon(genome_record)    # taxonomy id
    strain, crt = Strain.objects.get_or_create(taxonomy_id=strain_taxon,
                                               defaults={'name':org})
    # create genome object
    g = Genome(genome_accession=genome_record.name,
               #genome_type=None,
               sequence=genome_record.seq.tostring(),
               GC_content=gc,
               strain=strain)
    g.save()

def make_all_genes(genome_record):
    """Given genome record from NCBI db, create gene objects for all genes in
    the genome"""
    genes = bioutils.get_genes(genome_record)
    genome = Genome.objects.get(genome_accession=genome_record.name)
    assert genome
    for g in genes:
        Gene(genome=genome, **g).save()

def make_TF_instance(TF_rec):
    tf = TFInstance(protein_accession=TF_rec.name,
                    name=TF_rec.name,
                    description=TF_rec.description)
    tf.save()
    
