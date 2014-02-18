"""Template tags to render gene and regulation diagrams"""

from django.template import Library
from django.utils.safestring import mark_safe
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm

register = Library()

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

    feature = SeqFeature(FeatureLocation(site_instance.start, site_instance.end), strand=site_instance.strand)
    gds_features.add_feature(feature,
                             color=colors.red,
                             name="site",
                             label=False,
                             label_size=12)
    gdd.draw(format='linear',
             fragments=1,
             start=min(map(lambda r: r.gene.start, regulations)) - 150,
             end = max(map(lambda r: r.gene.end, regulations)) + 150,
             pagesize = (2*cm, 20*cm))
    return mark_safe(gdd.write_to_string('svg'))
