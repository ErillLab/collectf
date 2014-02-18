import re
from django import template
from django.utils.safestring import mark_safe

register = template.Library()

@register.filter
def HTMLify(description_text):
    """Given a description text, check if it has special text that needs to be
    processed for proper rendering."""
    # replace PMIDs
    description_text = re.sub("\[PMID::([0-9]+)\]",
                              '<sup><a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=\g<1>">[\g<1>]</a></sup>',
                              description_text)
    # replace PFAMs
    description_text = re.sub("\[PFAM::(\w+)\]",
                              '<sup><a href="http://pfam.sanger.ac.uk/family/\g<1>">[\g<1>]</a></sup>',
                              description_text)
    return mark_safe(description_text)
    
