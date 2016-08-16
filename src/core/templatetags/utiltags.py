"""
This file contains custom template tags for general use in the project.

See Django documentation on template tags:
https://docs.djangoproject.com/en/dev/ref/templates/builtins/
"""

import re

from django.utils.safestring import mark_safe
from django.template import Library
register = Library()

"""
    creates a new Django directive (for use in Jinja; Django's HTML)
    in this case, this allows getting an item from a dictionary, given the key
    e.g. <td>{{ motif_report.regulation_stats|get_item:"not specified" }}%</td> in view_motif_reports.html
         the format here is equivalent to get_item(regulation_stats, "not specified")
"""
@register.filter
def get_item(dictionary, key):
    """Get the dictionary element having the specified key. Lets use
    dictionaries inside templates."""
    return dictionary.get(key)


@register.filter
def get(o, index):
    return o[index]


@register.filter
def get_keys(dictionary):
    """Return all keys of a dictionary."""
    return sorted(dictionary.keys())


@register.filter
def HTMLify(description_text):
    """Turns the description text into HTML."""
    # replace PMIDs
    description_text = re.sub(
        '\[PMID::([0-9]+)\]',
        '<sup><a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=\g<1>">[\g<1>]</a></sup>',
        description_text)
    # replace PFAMs
    description_text = re.sub(
        '\[PFAM::(\w+)\]',
        '<sup><a href="http://pfam.sanger.ac.uk/family/\g<1>">[\g<1>]</a></sup>',
        description_text)
    return mark_safe(description_text)
