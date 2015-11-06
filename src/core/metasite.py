from . import models

class MetaSite:
    """MetaSite class definition.
    
    A site could be reported in multiple papers. To avoid redundancy
    (i.e. presenting the same site sequence multiple times, one per paper), the
    Curation_SiteInstance objects are collapsed into one entity called a
    "meta-site". To be collapsed into the same meta-site, two sites are
    required to overlap more than a defined threshold.
    """

    def __init__(self, curation_site_instances):
        self.curation_site_instances = curation_site_instances

