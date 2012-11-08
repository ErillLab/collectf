from django import template
from django.template import Context
from django.template.loader import get_template
from django.utils.safestring import mark_safe

register = template.Library()

@register.filter
def as_publication(pub):
    template = get_template("pub_template.html")
    c = Context({'pub': pub})
    return template.render(c)

