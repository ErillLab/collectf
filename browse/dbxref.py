"""Module for dbxref identifier conversion."""

def to_uniprot_dbxref(tf_instance_id):
    """Converts TFInstance id to UniProt dbxref identifier."""
    return 'EXPREG_' + '00' + hex(int(tf_instance_id))[2:].zfill(5) + '0'

def from_uniprot_dbxref(dbxref_id):
    """Converts UniProt dbxref identifier to CollecTF TFInstance id."""
    return int(dbxref_id[2:-1], 16)
