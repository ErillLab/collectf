"""Custom template tags for curation form."""

from django import template
from django.template import Context
from django.template.loader import get_template

register = template.Library()

@register.filter
def print_pub(pub):
    """Renders publication template."""
    templ = get_template("pub_template.html")
    c = Context({'pub': pub})
    return templ.render(c)

