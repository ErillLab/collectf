from Bio import Entrez

Entrez.email = "sefakilic@gmail.com"

def get_pubmed(pmid):
    """Retrieve pubmed publication from NCBI database."""
    handle = Entrez.esummary(db="pubmed", id=pmid)
    if handle:
        record = Entrez.read(handle)
        return record[0]
    
