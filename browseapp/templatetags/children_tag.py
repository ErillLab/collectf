from django.template import Library
register = Library()

@register.inclusion_tag('tax_children.html')
def children_tag(tax):
    children = tax.taxonomy_set.all()
    return {'children': children}

@register.filter
def has_children(tax):
    return bool(tax.taxonomy_set.all())
