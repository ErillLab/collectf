from django.template import Library

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

@register.filter
def my_add(sa, sb):
    # built-in add templatetag doesn't concat number and string
    return str(sa)+str(sb)

@register.filter
def defer_field(queryset, field):
    """Ignore field when retrieving object"""
    return queryset.defer(field)

@register.filter
def unique_TF_names(TF_list):
    # given a list of TF names, remove duplicates
    return list(set(lambda TF: TF.name, TF_list))
