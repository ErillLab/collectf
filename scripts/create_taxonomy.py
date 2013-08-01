from collectfapp import models
from collectfapp import bioutils
from Bio import Entrez
from Bio import SeqIO
import time
Entrez.email = 'sefakilic@gmail.com'

def create_taxonomy(genome):
    time.sleep(3)
    print 'creating taxonomy for', genome.genome_accession
    h = Entrez.efetch(db='nuccore', id=genome.genome_accession, retmode='gbwithparts', rettype='text')
    seq_record = SeqIO.read(h, 'gb')
    h.close()

    r = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='nuccore', id=seq_record.annotations['gi'], linkname='nuccore_taxonomy'))
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
        p,_ = models.Taxonomy.objects.get_or_create(rank=item['Rank'],
                                                  taxonomy_id=item['TaxId'],
                                                  name=item['ScientificName'],
                                                  parent=p)
    
    genome.taxonomy = p
    
    # add chromosome, organism and gi once we download genbank file
    gi = seq_record.annotations['gi']
    organism = seq_record.annotations['organism']
    
    assert seq_record.features[0].type=='source'
    q = seq_record.features[0].qualifiers
    chromosome = q['chromosome'][0] if 'chromosome' in q else ''

    genome.gi = gi
    genome.organism = organism
    genome.chromosome = chromosome
    genome.save()
    
def create_taxonomy_all():
    # For all strains in the genome, create all taxonomy levels up to Phylum.
    for genome in models.Genome.objects.iterator():
        create_taxonomy(genome)
        
def run():
    create_taxonomy_all()
    #genome = models.Genome.objects.get(genome_accession='NC_020286.1')
    #create_taxonomy(genome)
