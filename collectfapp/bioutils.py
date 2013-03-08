from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import time

Entrez.email = "sefakilic@gmail.com"

def get_pubmed(pmid):
    """Retrieve pubmed publication from NCBI database."""
    try:
        handle = Entrez.esummary(db="pubmed", id=pmid)
        record = Entrez.read(handle)
        return record[0]
    except RuntimeError:
        print 'err'
        return None

def get_seq_rec(accession, db):
    """Base function to retrieve genome or TF sequence record from NCBI database"""
    try:
        h = Entrez.efetch(db=db, rettype="gb", retmode="text", id=accession)
        seq_record = SeqIO.read(h, "gb")
        h.close()
        return seq_record
    except:
        print 'err2'
        return None

def get_genome(accession):
    """Retrieve genome from NCBI database"""
    return get_seq_rec(accession, db="nuccore")

def get_TF(accession):
    """Retrieve transcription factor from NCBI database"""
    return get_seq_rec(accession, db="protein")


def get_gene_id(feature):
    # extract gene id
    i = 0
    while not feature.qualifiers['db_xref'][i].startswith('GeneID:'):
        i += 1
    return feature.qualifiers['db_xref'][i][7:]
    
def get_genes(genome_rec):
    """Get list of all genes"""
    genes = [] # return list of genes
    # Use Entrez post method, because get method has limitation on url length
    # use Epost to post list of ids first
    # get gene ids
    gene_features = [f for f in genome_rec.features if f.type == 'gene']
    gids = [get_gene_id(f) for f in gene_features]
    # Occasionally, when collecTF tries to retrieve all gene summaries, NCBI refuses
    # to return all them. There should be some sort of limit for a query. Therefore,
    # the gene list summary query is chunked into 1000 gene pieces.
    recs = []
    chunk_size = 10
    for sliced_gids in [gids[i:i+chunk_size] for i in range(0, len(gids), chunk_size)]:
        print 'getting slice'
        req = Entrez.epost('gene', id=','.join(sliced_gids))
        res = Entrez.read(req)
        recs.extend(Entrez.read(Entrez.esummary(db='gene', webenv=res['WebEnv'],
                                                query_key=res['QueryKey'])))
        time.sleep(1)

    # two sources of gene data: Entrez epost and gene features from genome record
    for gid, feat, rec in zip(gids, gene_features, recs):
        assert rec['Id'] == gid
        genes.append({'gene_accession': gid,
                      'name': rec['Name'],
                      'description': rec['Description'],
                      'start': feat.location.start.position,
                      'end': feat.location.end.position,
                      'strand': feat.strand,
                      'locus_tag': ','.join(feat.qualifiers['locus_tag'])})
    return genes
   
def get_org_name(genome_record):
    """Given genome record from NCBI db, get organism name"""
    return genome_record.annotations['organism']

def get_org_taxon(genome_record):
    """Given genome rec from NCBI db, get organism taxonomy id"""
    org = get_org_name(genome_record)
    handle = Entrez.esearch(db='taxonomy', term=org)
    rec = Entrez.read(handle)
    assert int(rec['Count']) == 1
    return rec['IdList'][0]
