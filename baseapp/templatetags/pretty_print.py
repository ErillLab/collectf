from django.template import Library
from django.utils.safestring import mark_safe
import gene_diagram

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

def print_site_match(reported_site, m, is_exact):
    """Given a match object, make the html snippet to display it nice.  I am not sure
    this is a proper solution (HTML in python), couldn't find a better&easier way to
    do though.
    """
    strand = '+' if m.match.strand==1 else '-'
    nearby_genes = ['%s (%s)' % (g.locus_tag, g.name) for g in m.nearby_genes]
    
    s = ""
    if is_exact:
        s += ('<span class="sequence"> %s %s(%d, %d)</span><br/>' %
              (m.match.seq, strand,  m.match.start, m.match.end))
    else:
        s += print_alignment(reported_site, m.match.seq)

    s += (gene_diagram.site_match_diagram(m) +
          '<table class="table table-condensed">' +
          '<thead><tr><th>locus tag</th><th>gene name</th><th>function</th></tr></thead>' +
          '<tbody>'
          )
    for g in m.nearby_genes:
        s += "<tr><td>%s</td><td>%s</td><td>%s</td></tr>" % (g.locus_tag, g.name, g.description)
    s += ('</tbody>' + "</table><br/>")

    return mark_safe(s)  # render newline correctly
    

def print_alignment(seqa, seqb):
    """Given two sequences, pairwise align them and output HTML for curation
    exact/inexact site match steps"""
    assert len(seqa) == len(seqb)
    s = []
    s.append('<span class="sequence">')
    s.append("%s<br/>" % seqa)
    s.append(''.join('|' if seqa[i]==seqb[i] else '&nbsp;' for i in range(len(seqa))) + '<br/>')
    s.append('%s</br>' % seqb)
    s.append('</span>')
    return mark_safe( ''.join(s))
