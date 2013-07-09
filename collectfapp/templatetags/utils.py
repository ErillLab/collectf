from django.template import Library
from django.utils.safestring import mark_safe
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
import re
import operator

register = Library()

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter
def get_keys(dictionary):
    return dictionary.keys()


@register.filter
def get_keys_sorted_by_val(dictionary):
    sorted_dict = sorted(dictionary.iteritems(), key=operator.itemgetter(1))
    return map(lambda x: x[0], sorted_dict)
