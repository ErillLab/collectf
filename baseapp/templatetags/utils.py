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
