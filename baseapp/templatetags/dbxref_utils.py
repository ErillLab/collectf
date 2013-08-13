from django.template import Library

register = Library()

@register.filter
def id2dbxref(id):
    return 'CollecTF:EXPSITE_' + '00' + hex(int(id))[2:].zfill(5) + '0'

@register.filter
def id2dbxref_hex_only(id):
    return '00' + hex(int(id))[2:].zfill(5) + '0'
    

@register.filter
def dbxref2id(dbxref):
    return int(dbxref[2:-1], 16)
