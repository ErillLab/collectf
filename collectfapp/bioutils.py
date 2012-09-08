from Bio import Entrez

Entrez.email = "sefakilic@gmail.com"

def get_pubmed(pmid):
    """Retrieve pubmed publication from NCBI database."""
    try:
        handle = Entrez.esummary(db="pubmed", id=pmid)
        record = Entrez.read(handle)
        return record[0]
    except RuntimeError:
        return None
    
