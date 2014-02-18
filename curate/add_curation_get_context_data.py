"""For each step of the curation submission, return the step-specific context data."""
import session_utils

def publication_context_data(wiz):
    return {}

def genome_context_data(wiz):
    return {}

def techniques_context_data(wiz):
    return {}

def site_entry_context_data(wiz):
    return {}

def site_exact_match_context_data(wiz):
    return {}

def site_soft_match_context_data(wiz):
    return {}

def site_annotation_context_data(wiz):
    sites = session_utils.get(wiz.request.session, "sites")
    return {"sites": sites}
        
def gene_regulation_context_data(wiz):
    return {}

def curation_review_context_data(wiz):
    return {}
