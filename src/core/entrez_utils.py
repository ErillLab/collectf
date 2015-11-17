"""Functions to fetch information from NCBI via Entrez."""

import urllib
import uniprot

from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'sefa1@umbc.edu'


class EntrezException(Exception):
    pass


def get_pubmed(pmid):
    """Retrieves the PubMed publication from NCBI database."""
    try:
        handle = Entrez.esummary(db='pubmed', id=pmid)
        record = Entrez.read(handle)
        return record[0]
    except RuntimeError:
        raise EntrezException


def get_genome(accession):
    """Retrieves the genome from RefSeq database."""
    try:
        handle = Entrez.efetch(db='nuccore', id=accession,
                               retmode='gbwithparts', rettype='text')
        record = SeqIO.read(handle, 'gb')
        handle.close()
        return record
    except urllib.error.HTTPError:
        raise EntrezException


def get_organism_taxon(genome_record):
    """Finds organism taxonomy ID of a given record using Elink utility."""
    try:
        gi = genome_record.annotations['gi']
        record = Entrez.read(Entrez.elink(
            db='taxonomy', dbfrom='nuccore', id=gi, linkname='nuccore_taxonomy'))
        tax_id = record[0]['LinkSetDb'][0]['Link'][0]['Id']
        return tax_id
    except:
        raise EntrezException

def get_taxonomy(genome_record):
    """Retrieves the taxonomy record from NCBI."""
    link_record = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='nuccore',
                                           id=genome_record.annotations['gi'],
                                           linkname='nuccore_taxonomy'))
    tax_id = link_record[0]['LinkSetDb'][0]['Link'][0]['Id']
    handle = Entrez.efetch(db='Taxonomy', id=tax_id, retmode='xml')
    record = Entrez.read(handle)[0]
    # Add the taxon itself to LineageEx.
    record['LineageEx'].append({'Rank': record['Rank'],
                                'ScientificName': record['ScientificName'],
                                'TaxId': record['TaxId']})
    
    return record



def get_genes(genome_record):
    """Gets the list of all genes, given a genome record object."""

    def get_gene_annotation(id_list):
        """Gets gene annotations.
        
        Uses Bio.Entrez epost to submit the data to NCBI, and esummary to
        retrieve the information. Returns a list gene summary objects.
        """
        epost_result = Entrez.read(Entrez.epost('gene', id=','.join(id_list)))

        # NCBI may refuse to return all gene summaries due to some sort of
        # query limit. Query summaries for 1000 genes at a time.
        
        runtime_error = 0  # error couunt  during Entrez esummary
        while runtime_error < 10: 
            try:
                request = Entrez.esummary(db="gene",
                                          webenv=epost_result["WebEnv"],
                                          query_key=epost_result["QueryKey"])
                records = Entrez.read(request)
            except RuntimeError,  e:
                print "Error occurred during epost+esummary:", e
                print "Trying again."
                runtime_error += 1
            else:
                break
        return records['DocumentSummarySet']['DocumentSummary']

    def get_gene_id(feature):
        """Extracts the gene id, given a Biopython SeqFeature object."""
        i = 0
        while not feature.qualifiers['db_xref'][i].startswith('GeneID:'):
            i += 1
        gene_id = feature.qualifiers['db_xref'][i][7:]
        return gene_id

    genes = [] # list of genes
    gene_features = [f for f in genome_record.features if f.type == 'gene']
    gids = [get_gene_id(f) for f in gene_features]
    recs = []
    chunk_size = 1000
    for start in xrange(0, len(gids), chunk_size):
        end = min(len(gids), start+chunk_size)
        recs = recs + get_gene_annotation(gids[start:end])

    # two sources of gene data: Entrez epost and gene features from genome
    # record
    for gid, feat, rec in zip(gids, gene_features, recs):
        assert rec.attributes['uid'] == gid
        genes.append({'gene_accession': gid,
                      'name': rec['Name'],
                      'description': rec['Description'],
                      'start': feat.location.start.position,
                      'end': feat.location.end.position,
                      'strand': feat.strand,
                      'locus_tag': ','.join(feat.qualifiers['locus_tag'])})
    return genes


def get_uniprot_TF(accession):
    """Retrieves UniProt record for the given accession number."""
    return uniprot.retrieve(accession)


def get_refseq_TF(accession):
    """Retrieve transcription factor from NCBI database."""
    try:
        handle = Entrez.efetch(
            db='protein', id=accession, retmode='text', rettype='gb')
        seq_record = SeqIO.read(handle, 'gb')
        return seq_record
    except Exception:
        raise EntrezException
