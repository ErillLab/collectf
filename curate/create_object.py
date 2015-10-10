"""List of functions to create objects in the database."""

from Bio import Entrez
from Bio import SeqIO

from base import bioutils
from models import *

def make_pub(pubrec, cd):
    """Creates models.Publication object given the publication record."""
    pmid = cd.get("pmid", None)  # None if not pubmed publication
    url ="http://www.ncbi.nlm.nih.gov/pubmed?term=%s" % pmid if pmid else cd['URL']
    publication_type = "pubmed" if pmid else "nonpubmed"
    return Publication(publication_type=publication_type,
                       pmid=pmid,
                       authors=', '.join(pubrec.get("AuthorList")),
                       journal=pubrec.get("FullJournalName"),
                       title=unicode(pubrec.get("Title")),
                       publication_date=pubrec.get("PubDate"),
                       volume=pubrec.get("Volume"),
                       issue=pubrec.get("Issue"),
                       pages=pubrec.get("Pages"),
                       url=url,
                       pdf=None,
                       reported_TF=cd["reported_TF"],
                       reported_species=cd["reported_species"],
                       contains_promoter_data=cd["contains_promoter_data"],
                       contains_expression_data=cd["contains_expression_data"],
                       submission_notes=cd["submission_notes"],
                       curation_complete=False)

def create_taxonomy(genome_record):
    """Creates all taxonomy items up to Bacteria given the genome record.

    Returns the most specific rank.
    """
    # Get taxonomy object from NCBI
    rec = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='nuccore',
                                   id=genome_record.annotations['gi'],
                                   linkname='nuccore_taxonomy'))
    tax_id = rec[0]['LinkSetDb'][0]['Link'][0]['Id']
    handle = Entrez.efetch(db='Taxonomy', id=tax_id, retmode='xml')
    records = Entrez.read(handle)
    assert len(records) == 1, "More than one taxonomy records retrieved."
    record = records[0]
    record['LineageEx'].append({'Rank': records[0]['Rank'],
                                'ScientificName': records[0]['ScientificName'],
                                'TaxId': records[0]['TaxId']})

    # The lineage that NCBI returns seems already sorted, but to be safe, sort
    # it here.
    lineage = []
    for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
        r = filter(lambda x: x['Rank']==rank, record['LineageEx'])
        assert len(r) <= 1
        lineage.extend(r)
    assert len(lineage) >= 1, "Number of lineages greater than one"
    p = None
    for item in lineage:
        p,_ = Taxonomy.objects.get_or_create(
            rank=item['Rank'], taxonomy_id=item['TaxId'],
            name=item['ScientificName'], parent=p)
    return p

def make_genome(genome_record, strain_tax):
    """Creates a models.Genome object given the genome record."""
    gc = bioutils.GC(genome_record.seq)  # GC content
    p = create_taxonomy(genome_record)   # Create all taxonomy
    assert genome_record.features[0].type == 'source'
    q = genome_record.features[0].qualifiers
    chromosome = q['chromosome'][0] if 'chromosome' in q else ''
    # Create genome sequence object
    genome_sequence = GenomeSequence.objects.create(
        sequence=genome_record.seq.tostring())
    # Create genome object
    genome = Genome(genome_accession=genome_record.id,
                    genome_sequence=genome_sequence,
                    GC_content=gc,
                    gi=genome_record.annotations['gi'],
                    organism=genome_record.annotations['organism'],
                    chromosome=chromosome,
                    taxonomy=p)

    genes = bioutils.get_genes(genome_record)
    if genes:
        genome.save()
        _ = [Gene(genome=genome, **gene).save() for gene in genes]
        return genome
    return None

def make_all_genes(genome_record, genome_obj):
    """Creates gene objects for all genes in the genome."""
    genes = bioutils.get_genes(genome_record)
    return [Gene(genome=genome_obj, **g) for g in genes]

def make_TF_instance(TF_rec, TF):
    """Creates the TF object in the database, given the TF record."""
    TF_instance = TFInstance(
        protein_accession=TF_rec.name, name=TF_rec.name,
        description=TF_rec.description, TF=TF)
    TF_instance.save()
