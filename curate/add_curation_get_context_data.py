"""For each step of the curation submission, return the step-specific context
data."""
import curate.session_utils as session_utils
from base import bioutils

def publication_context_data(wiz):
    return {}

def genome_context_data(wiz):
    return {}

def techniques_context_data(wiz):
    return {}

def site_entry_context_data(wiz):
    return {}

def site_exact_match_context_data(wiz):
    # Generate the weblogo of the entered sites
    sites = session_utils.get(wiz.request.session, 'sites')
    site_type = session_utils.get(wiz.request.session, 'site_type')
    d = {}
    if (site_type not in ['non_motif_associated', 'var_motif_associated'] and (len(sites)>1)):
        d['weblogo_img'] = bioutils.weblogo_uri([site.seq for site in sites])
    return d

def site_soft_match_context_data(wiz):
    d = {}
    # If no soft match, put something to the context data, so that client-side
    # knows that there is nothing in this page and skip it accordingly.
    sites = session_utils.get(wiz.request.session, 'sites')
    if not any(site.is_matched() and not site.get_match().is_exact()
               for site in sites):
        d['no_soft_match'] = True
    return d

def site_annotation_context_data(wiz):
    sites = session_utils.get(wiz.request.session, "sites")
    techniques = session_utils.get(wiz.request.session, "techniques")
    has_qdat = session_utils.get(wiz.request.session, 'has_quantitative_data')
    return {"sites": sites, #[site for site in sites if site.is_matched()],
            "techniques": techniques,
            "has_quantitative_data": has_qdat,
            }
        
def gene_regulation_context_data(wiz):
    return {}

def curation_review_context_data(wiz):
    return {}
