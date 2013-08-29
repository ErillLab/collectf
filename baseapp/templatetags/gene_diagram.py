from django.template import Library
from django.utils.safestring import mark_safe
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm

register = Library()

@register.filter
def site_match_diagram(site_match):
    """This method is called during curation submission. When the reported sites are
    given, they are searched in the genome and matched are displayed. For display,
    site match diagram is created and presented."""
    gdd = GenomeDiagram.Diagram("Site match diagram")
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()
    # draw genes
    for gene in site_match.nearby_genes:
        feature = SeqFeature(FeatureLocation(gene.start, gene.end), strand=gene.strand)
        gds_features.add_feature(feature,
                                 name=gene.name,
                                 label=True,
                                 label_size=10,
                                 label_angle= 0 if gene.strand==1 else 180,
                                 label_position='middle',
                                 sigil='ARROW',
                                 arrowshaft_height=1.0,
                                 color=colors.lightblue)
    # draw binding site
    feature = SeqFeature(FeatureLocation(site_match.match.start,site_match.match.end),
                         strand=site_match.match.strand)
    gds_features.add_feature(feature,
                             color=colors.red,
                             name='site',
                             label=False,
                             label_size=12)
    gdd.draw(format='linear',
             fragments=1,
             start=min(map(lambda g: g.start, site_match.nearby_genes))-50,
             end=max(map(lambda g: g.end, site_match.nearby_genes))+50,
             pagesize=(3*cm, 20*cm))
    return mark_safe(gdd.write_to_string('svg'))
                         
@register.filter
def regulation_diagram(regulations, site_instance):
    """This method is called for drawing the regulation diagram."""
    gdd = GenomeDiagram.Diagram('Regulation Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()
    # draw genes
    for reg in regulations:
        feature = SeqFeature(FeatureLocation(reg.gene.start, reg.gene.end), strand=reg.gene.strand)
        gds_features.add_feature(feature,
                                 name=reg.gene.name,
                                 label=True,
                                 label_size=10,
                                 label_angle=0 if reg.gene.strand==1 else 180,
                                 label_position='middle',
                                 sigil="ARROW",
                                 arrowshaft_height=1.0,
                                 color=colors.green if reg.evidence_type=="exp_verified" else colors.grey)
    # draw binding site

    # site_instance argument may be type of SiteInstance or MetaSiteInstance. If
    # MetaSiteInstance, there is no strand, set to 1.
    if not hasattr(site_instance, 'strand'):
        site_instance.strand = 1

    feature = SeqFeature(FeatureLocation(site_instance.start, site_instance.end), strand=site_instance.strand)
    gds_features.add_feature(feature,
                             color=colors.red,
                             name="site",
                             label=False,
                             label_size=12)
    gdd.draw(format='linear',
             fragments=1,
             start=min(map(lambda r: r.gene.start, regulations)) - 50,
             end = max(map(lambda r: r.gene.end, regulations)) + 50,
             pagesize = (2*cm, 20*cm))
    return mark_safe(gdd.write_to_string('svg'))
