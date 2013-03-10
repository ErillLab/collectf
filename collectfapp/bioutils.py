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
    gene_id = feature.qualifiers['db_xref'][i][7:]
    return gene_id

def get_gene_annotation(id_list):
    """Annotates Entrez Gene ids using Bio.Entrez, in particular epost (to submit the
    data to NCBI) and esummary to retrieve the information. Returns a list gene
    summary objects."""
    epost_result = Entrez.read(Entrez.epost("gene", id=','.join(id_list)))
    # Occasionally, when collecTF tries to retrieve all gene summaries, NCBI refuses
    # to return all them. There should be some sort of limit for a query. Therefore,
    # the gene list summary query is chunked into 1000 gene pieces.
    
    # edit: After several different ways to fix it, the best way is breaking up the
    # list into batches of size 1000 AND retry if any batch fails for some reason.
    runtime_error = 0  # number of runtime errors during Entrez esummary
    
    while runtime_error < 10:  # if runtime error is consistent, there is no point trying again
        try:
            request = Entrez.esummary(db="gene", webenv=epost_result["WebEnv"],
                                      query_key=epost_result["QueryKey"])
            records = Entrez.read(request)
        except RuntimeError as e:
            print "Error occurred during epost+esummary:", e
            print "Trying again."
            runtime_error += 1
        else:
            runtime_error = 0
            break

    assert runtime_error == 0
    return records


def get_genes(genome_rec):
    """Get list of all genes"""
    genes = [] # return list of genes
    # Use Entrez post method, because get method has limitation on url length
    # use Epost to post list of ids first
    # get gene ids
    gene_features = [f for f in genome_rec.features if f.type == 'gene']
    gids = [get_gene_id(f) for f in gene_features]
    
    recs = []
    chunk_size = 1000
    for start in xrange(0, len(gids), chunk_size):
        end = min(len(gids), start+chunk_size)
        #print gids[start:end]
        recs = recs + get_gene_annotation(gids[start:end])
        
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
