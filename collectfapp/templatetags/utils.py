from django.template import Library
from django.utils.safestring import mark_safe
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm
import re


register = Library()

@register.filter
def get_item(dictionary, key):
    return dictionary.get(key)

@register.filter
def get_keys(dictionary):
    return dictionary.keys()

@register.filter
def regulation_diagram(regulations, site_instance):
    """Given a site instance and regulation objects, output HTML diagram"""
    gdd = GenomeDiagram.Diagram('Regulation Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    # draw genes
    for reg in regulations:
        feature = SeqFeature(FeatureLocation(reg.gene.start, reg.gene.end), strand=reg.gene.strand)
        label_angle = 0 if reg.gene.strand == 1 else 180
        label_position = 'middle'
        gds_features.add_feature(feature, name=reg.gene.name, label=True,
                                 label_size=15, label_angle=label_angle,
                                 label_position=label_position,
                                 sigil="ARROW", arrowshaft_height=1.0,
                                 color=colors.green if reg.evidence_type=="exp_verified" else colors.grey)
        
    # draw binding site
    feature = SeqFeature(FeatureLocation(site_instance.start, site_instance.end), strand=site_instance.strand)
    gds_features.add_feature(feature, color=colors.red, name="site", label=False, label_size=12)

    gdd.draw(format='linear', fragments=1,
             start=min(map(lambda r: r.gene.start, regulations)) - 50,
             end = max(map(lambda r: r.gene.end, regulations)) + 50,
             pagesize = (3*cm, 17*cm))
    
    return mark_safe(gdd.write_to_string('svg'))
    
def site_match_diagram(site_match):
    gdd = GenomeDiagram.Diagram('Site Match diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    # draw genes
    for g in site_match.nearby_genes:
        feature = SeqFeature(FeatureLocation(g.start, g.end), strand=g.strand)
        label_angle = 0 if g.strand == 1 else 180
        gds_features.add_feature(feature, name=g.name, label=True,
                                 label_size=14, label_angle=label_angle,
                                 label_position = 'middle',
                                 sigil="ARROW", arrowshaft_height=1.0,
                                 color=colors.lightblue)
    # draw binding site
    feature = SeqFeature(FeatureLocation(site_match.match.start, site_match.match.end),
                         strand=site_match.match.strand)
    gds_features.add_feature(feature, color=colors.red, name="site", label=True, label_size=12)
        
    gdd.draw(format='linear', fragments=1,
             start=min(map(lambda g: g.start, site_match.nearby_genes))-50,
             end = max(map(lambda g: g.end, site_match.nearby_genes))+50,
             pagesize = (3*cm, 17*cm))
    return mark_safe(gdd.write_to_string('svg'))


@register.filter
def description_process(description_text):
    """Given a description text, check if it has special text that needs to be
    processed."""

    # replace PMIDs
    description_text = re.sub("\[PMID::([0-9]+)\]",
                              '<sup><a href="http://www.ncbi.nlm.nih.gov/pubmed/?term=\g<1>">[\g<1>]</a></sup>',
                              description_text)
    # replace PFAMs
    description_text = re.sub("\[PFAM::(\w+)\]",
                              '<sup><a href="http://pfam.sanger.ac.uk/family/\g<1>">[\g<1>]</a></sup>',
                              description_text)
    return mark_safe(description_text)
    

