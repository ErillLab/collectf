from django.template import Library
from django.utils.safestring import mark_safe


register = Library()

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter
def get_keys(dictionary):
    return dictionary.keys()

@register.filter
def gene_diagram(regulations, site_instance):

    """Given a site instance and regulation objects, sort them and output HTML for the diagram"""
    regulations.sort(key=lambda r: r.gene.start)
    left_regulations = [r for r in regulations if r.gene.start < site_instance.start]
    right_regulations = [r for r in regulations if r.gene.start > site_instance.start]
    site_html_output = ('<div class="tf-img" data-toggle="popover" data-content="%s"></div>' %
                        site_instance.seq)
    return (gene_diagram_helper(left_regulations) +
            mark_safe(site_html_output) +
            gene_diagram_helper(right_regulations))
    
def gene_diagram_helper(regulations):
    """Draw regulated genes"""
    html_output = ""
    for reg in regulations:
        assert reg.evidence_type in ["exp_verified", "inferred"]
        assert reg.gene.strand in [1, -1]
        # check css for gene-img class
        gene_img_div = "gene-img %s %s" % ("pos" if reg.gene.strand == 1 else "neg",
                                           "exp" if reg.evidence_type == "exp_verified" else "inf")
        html_output += ('<div class="%s" data-toggle="popover" data-original-title="%s" data-content="%s"></div>' %
                        (gene_img_div, reg.gene.name, reg.gene.description))
        
    return mark_safe(html_output)
    
    
