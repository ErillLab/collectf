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
def regulation_diagram(regulations, site_instance):
    """Given a site instance and regulation objects, sort them and output HTML for the diagram"""
    def regulation_diagram_helper(regulations):
        """Draw regulated genes"""
        html_output = ""
        for reg in regulations:
            assert reg.evidence_type in ["exp_verified", "inferred"]
            assert reg.gene.strand in [1, -1]
            # check css for gene-img class
            gene_img_div = "gene-img %s %s" % ("pos" if reg.gene.strand == 1 else "neg",
                                               "exp" if reg.evidence_type == "exp_verified" else "inf")
            data_title = reg.gene.name
            data_content = gene_info(reg.gene)
            html_output += ('<div class="%s" data-toggle="popover" data-original-title="%s" data-content="%s"></div>' %
                            (gene_img_div, data_title, data_content))
        return mark_safe(html_output)

    regulations.sort(key=lambda r: r.gene.start)
    left_regulations = [r for r in regulations if r.gene.start < site_instance.start]
    right_regulations = [r for r in regulations if r.gene.start >= site_instance.start]
    site_html_output = ('<div class="tf-img" data-toggle="popover" data-content="%s"></div>' %
                        mark_safe(site_info(site_instance)))
    return (regulation_diagram_helper(left_regulations) +
            mark_safe(site_html_output) +
            regulation_diagram_helper(right_regulations))


def biopython_diagram(site_match):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Graphics import GenomeDiagram
    from reportlab.lib import colors
    from reportlab.lib.units import cm

    gdd = GenomeDiagram.Diagram('Test Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    for g in site_match.nearby_genes:
        feature = SeqFeature(FeatureLocation(g.start, g.end), strand=g.strand)
        gds_features.add_feature(feature, name=g.name, label=True,
                                 sigil="ARROW", arrowshaft_height=1.0,
                                 color=colors.lightblue)
        
    gdd.draw(format='linear', fragments=1,
             start=min(map(lambda g: g.start, site_match.nearby_genes))-50,
             end = max(map(lambda g: g.end, site_match.nearby_genes))+50)


    return gdd.write_to_string('svg')


def match_diagram(site_match):
    """Ouput HTML for the diagram of site search matches.  The difference from the
    previous one is that this method is called on curation forms. Since there is no
    regulation information at that time, all genes are visualized as _gray_ block
    arrows."""
    def gene_diagram_helper(genes):
        """Draw a set of genes"""
        html_output = ""
        for g in genes:
            assert g.strand in [1, -1]
            gene_img_div = "gene-img match %s" % ("pos" if g.strand == 1 else "neg")
            data_title = g.name
            data_content = gene_info(g)
            html_output += ('<div class="%s" data-toggle="popover" data-original-title="%s" data-content="%s"></div>' %
                            (gene_img_div, data_title, data_content))
        return mark_safe(html_output)

    site_match.nearby_genes.sort(key=lambda gene: gene.start)
    left_genes = [g for g in site_match.nearby_genes if g.start < site_match.match.start]
    right_genes = [g for g in site_match.nearby_genes if g.start >= site_match.match.start]
    data_content = site_info(site_match.match)
    site_html_output = ('<div class="tf-img" data-toggle="popover" data-content="%s"></div>' %
                        mark_safe(data_content))
    html_out = mark_safe('<br/>' + gene_diagram_helper(left_genes) +
                         mark_safe(site_html_output) +
                         gene_diagram_helper(right_genes) + '<br/>')
    return html_out

def site_info(site_instance):
    """Given a site instance, produce info in HTML"""
    data_content = '''<b>site:</b> %s<br/>
                      <b>location: </b>%s (%d, %d)</br>''' % (
                      site_instance.seq,
                      '+' if site_instance.strand == 1 else '-',
                      site_instance.start, site_instance.end)
    return data_content
                      
def gene_info(gene):
    """Given a gene, produce information for that gene in HTML"""
    data_content = ("""<b>name:</b> %s<br/>
                       <b>locus tag:</b> %s<br/>
                       <b>location:</b> %s (%d, %d)<br/>
                       <b>function:</b> %s""" %
                    (gene.name, gene.locus_tag,
                     '+' if gene.strand == 1 else '-',
                     gene.start, gene.end, gene.description))

    return data_content

    

