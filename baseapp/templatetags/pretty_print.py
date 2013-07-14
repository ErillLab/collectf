from django.template import Library
from django.utils.safestring import mark_safe

register = Library()

@register.filter
def site2label(sid,site):
    return mark_safe('<span class="sequence">' +'[%d] %s' % (sid,site) + '</span>')

@register.filter
def match2label(id,match):
    if (len(match.match.seq) <= 25):
        label = match.match.seq
        
    else:
        strand = '+' if match.match.strand==1 else '-'
        label = 'loc %s[%d,%d]' % (strand, match.match.start, match.match.end)

    label =  ('[%d] %s' % (id,label))
    return mark_safe(label)

