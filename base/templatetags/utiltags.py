"""
This file contains custom template tags for general use in the project.

See Django documentation on template tags:
https://docs.djangoproject.com/en/dev/ref/templates/builtins/
"""

from django.template import Library
register = Library()

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

