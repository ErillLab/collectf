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
    return sorted(dictionary.keys())


@register.filter
def sort_by(queryset, order):
    return queryset.order_by(order)
