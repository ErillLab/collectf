"""Functions to fetch information from NCBI via Entrez."""

import uniprot

from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'erill@umbc.edu'


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
                               rettype='gbwithparts', retmode='text')
        record = SeqIO.read(handle, 'gb')
        handle.close()
        return record
    except:
        raise EntrezException


def get_organism_taxon(genome_record):
    """Finds organism taxonomy ID of a given record using Elink utility."""
    try:
        src_feature = genome_record.features[0]
        db_xrefs = src_feature.qualifiers['db_xref']
        taxon_db_xref, = [db_xref for db_xref in db_xrefs
                          if db_xref.startswith('taxon:')]
        tax_id = taxon_db_xref.replace('taxon:', '')
        print tax_id
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
    genes = []                  # list of genes
    feature_index = 0
    features = genome_record.features
    while feature_index < len(features):
        if features[feature_index].type == 'gene':
            gene_feature = features[feature_index]
            gene_rec = features[feature_index].qualifiers
            locus_tag = gene_rec['locus_tag']
            # Check if there is a product (CDS, tRNA, etc. ) entry for the gene
            # and get the description, if possible.
            gene_type = ""
            description = ""
            protein_id = None
            if feature_index+1 < len(features):
                next_feature = features[feature_index+1]
                next_rec = features[feature_index+1].qualifiers
                if next_rec.get('locus_tag') == locus_tag:
                    description = ', '.join(next_rec.get('product', []))
                    gene_type = next_feature.type
                    protein_id = next_rec.get('protein_id')
                    feature_index += 1
                else:
                    print "No product for", gene_rec

            genes.append({'name': ', '.join(gene_rec.get('gene', locus_tag)),
                          'description': description,
                          'start': gene_feature.location.start.position,
                          'end': gene_feature.location.end.position,
                          'strand': gene_feature.strand,
                          'locus_tag': ', '.join(locus_tag),
                          'gene_type': gene_type,
                          # Additional stuff
                          'protein_id': protein_id,
                          'old_locus_tag': gene_rec.get('old_locus_tag', [])})

        feature_index += 1
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
