"""
This module contains the class implementation for motif reports.
"""

from django.utils.functional import cached_property

from core import bioutils
from core import metasite
from core import lasagna


class MotifReport:
    """Collection of meta-sites that have the same TF and species."""
    def __init__(self, curation_site_instances):
        self.meta_sites = metasite.create_meta_sites(curation_site_instances)

    def id(self):
        """Returns the unique identifier for the instance."""
        return id(self)

    @property
    def TF(self):
        """Returns the TF."""
        return self.meta_sites[0].delegate.curation.TF

    @property
    def TF_name(self):
        """Returns the name of the TF."""
        return self.TF.name.strip()

    @cached_property
    def TF_accessions(self):
        """Returns the TF accession numbers."""
        return self.meta_sites[0].delegate.curation.TF_instances.values_list(
            'uniprot_accession', flat=True).distinct()

    @property
    def species_name(self):
        """Returns the name of the species."""
        return self.meta_sites[0].delegate.site_instance.genome.organism

    @cached_property
    def aligned_sites(self):
        """Aligns all binding sites using LASAGNA."""
        site_instances = [meta_site.delegate.site_instance
                          for meta_site in self.meta_sites
                          if meta_site.site_type == 'motif_associated']
        return lasagna.lasagna(site_instances)

    @property
    def weblogo(self):
        """Returns the weblogo URI for its sites."""
        return bioutils.weblogo_uri(self.aligned_sites)



class EnsembleMotifReport:
    """Collection of meta-sites that are requested to be viewed together."""
    def __init__(self, curation_site_instances):
        self.meta_sites = metasite.create_meta_sites(curation_site_instances)
