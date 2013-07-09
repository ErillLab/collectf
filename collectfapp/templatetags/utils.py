from django.template import Library
from django.utils.safestring import mark_safe
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
import re


register = Library()

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter
def get_keys(dictionary):
    return dictionary.keys()




