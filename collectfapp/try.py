from Bio import Entrez
from Bio import SeqIO
Entrez.email = 'sefakilic@gmail.com'


def get_genome(accession):
    """Retrieve genome from NCBI database"""
    try:
        print accession
        h = Entrez.efetch(db='nuccore', id=accession, retmode='gbwithparts', rettype='text')
        seq_record = SeqIO.read(h, 'gb')
        print seq_record
        h.close()
        return seq_record
    except:
        print 'Something went wrong during refseq retrieval'
        return None
