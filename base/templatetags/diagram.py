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
def regulation_diagram(regulations):
    """Draws the regulation diagram of the given regulation."""
    gdd = GenomeDiagram.Diagram('Regulation Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()
    # Draw genes.
    for regulation in regulations:
        feature = SeqFeature(
            FeatureLocation(regulation.gene.start, regulation.gene.end),
            strand=regulation.gene.strand)
        gds_features.add_feature(
            feature,
            name=regulation.gene.name,
            label=True,
            label_size=10,
            label_angle=0 if regulation.gene.strand == 1 else 180,
            label_position='middle',
            sigil="ARROW",
            arrowshaft_height=1.0,
            color=(colors.green if regulation.evidence_type == 'exp_verified'
                   else colors.grey))

    # Draw binding site
    site_instance = regulations[0].curation_site_instance.site_instance
    feature = SeqFeature(
        FeatureLocation(site_instance.start, site_instance.end),
        strand=site_instance.strand)
    gds_features.add_feature(feature,
                             color=colors.red,
                             name='site',
                             label=False,
                             label_size=12)
    gdd.draw(format='linear',
             fragments=1,
             start=min(map(lambda r: r.gene.start, regulations)) - 150,
             end=max(map(lambda r: r.gene.end, regulations)) + 150,
             pagesize = (2*cm, 12*cm))
    return mark_safe(gdd.write_to_string('svg'))
