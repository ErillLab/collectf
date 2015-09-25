# Find RefSeq accession to UniProt ID for all proteins in CollecTF.

import uniprot

from base import models

def get_all_proteins():
    """Returns the list of protein accessions in CollecTF."""
    return [protein.protein_accession
            for protein in models.TFInstance.objects.all()]

def refseq_to_uniprot(refseq_accessions):
    """Maps all RefSeq accessions to UniProt IDs."""
    return uniprot.map(refseq_accessions, f='P_REFSEQ_AC', t='ACC')

def run():
    refseq_accs = get_all_proteins()
    uniprot_mapping = refseq_to_uniprot(refseq_accs)
    for refseq_acc in refseq_accs:
        uniprot_accs = uniprot_mapping.get(refseq_acc, set())
        print refseq_acc, list(uniprot_accs)

    
