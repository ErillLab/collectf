# List of functions to create objects in db

from models import *
import bioutils
from Bio import Entrez
from Bio import SeqIO

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
                    authors = ', '.join(pubrec.get("AuthorList")), 
                    journal = pubrec.get("FullJournalName"),
                    title = pubrec.get("Title"),
                    publication_date = pubrec.get("PubDate"),
                    volume = pubrec.get("Volume"),
                    issue = pubrec.get("Issue"),
                    pages = pubrec.get("Pages"),
                    url=url,
                    pdf=None,
                    reported_TF=cd["reported_TF"],
                    reported_species=cd["reported_species"],
                    contains_promoter_data=cd["contains_promoter_data"],
                    contains_expression_data=cd["contains_expression_data"],
                    submission_notes=cd["submission_notes"],
                    curation_complete=False)
    return p

def create_taxonomy(genome_record):
    """Given genome record, create all taxonomy elements up to Bacteria.
    Return the most specific rank.
    """
    # get taxonomy object from NCBI
    r = Entrez.read(Entrez.elink(db='taxonomy',
                                 dbfrom='nuccore',
                                 id=genome_record.annotations['gi'],
                                 linkname='nuccore_taxonomy'))
    
    tax_id = r[0]['LinkSetDb'][0]['Link'][0]['Id']
    h = Entrez.efetch(db='Taxonomy', id=tax_id, retmode='xml')
    records = Entrez.read(h)
    assert len(records) == 1
    record = records[0]
    record['LineageEx'].append({'Rank': records[0]['Rank'],
                                'ScientificName': records[0]['ScientificName'],
                                'TaxId': records[0]['TaxId']})
    # The lineage that NCBI returns seems already sorted, but to be safe, sort it here.
    lineage = []
    for rank in ['phylum', 'class', 'order', 'family', 'genus', 'species']:
        r = filter(lambda x: x['Rank']==rank, record['LineageEx'])
        assert len(r) <= 1
        lineage.extend(r)
    print lineage
    assert len(lineage) >= 1
    p = None
    for item in lineage:
        p,_ = Taxonomy.objects.get_or_create(rank=item['Rank'],
                                             taxonomy_id=item['TaxId'],
                                             name=item['ScientificName'],
                                             parent=p)
    return p

def make_genome(genome_record, strain_taxon):
    """Given genome record that is fetched from NCBI, create a models.Genome object"""
    gc = bioutils.GC(genome_record.seq)  # GC content
    p = create_taxonomy(genome_record) # create all taxonomy
    # chromosome info?
    assert genome_record.features[0].type=='source'
    q = genome_record.features[0].qualifiers
    chromosome = q['chromosome'][0] if 'chromosome' in q else ''
    # create genome object
    g = Genome(genome_accession=genome_record.id,
               sequence=genome_record.seq.tostring(),
               GC_content=gc,
               gi=genome_record.annotations['gi'],
               organism=genome_record.annotations['organism'],
               chromosome=chromosome,
               taxonomy=p)
    g.save()

def make_all_genes(genome_record):
    """Given genome record from NCBI db, create gene objects for all genes in
    the genome"""
    genes = bioutils.get_genes(genome_record)
    genome = Genome.objects.get(genome_accession=genome_record.id)
    assert genome
    for g in genes:
        Gene(genome=genome, **g).save()

def make_TF_instance(TF_rec):
    tf = TFInstance(protein_accession=TF_rec.name,
                    name=TF_rec.name,
                    description=TF_rec.description)
    tf.save()
    
