"""
This module contains the class implementation for motif reports.
"""

from django.utils.functional import cached_property

from core import bioutils
from core import metasite
from core import lasagna


class MotifReportBase:
    """Base class for MotifReport and EnsembleMotifReport."""
    def __init__(self, curation_site_instances):
        self.meta_sites = metasite.create_meta_sites(curation_site_instances)

    def id(self):
        """Returns the unique identifier for the instance."""
        return id(self)

    @cached_property
    def aligned_sites(self):
        """Aligns all binding sites using LASAGNA."""
        site_instances = [meta_site.delegate.site_instance
                          for meta_site in self.meta_sites
                          if meta_site.site_type == 'motif_associated']
        return lasagna.lasagna(site_instances)

    @cached_property
    def weblogo(self):
        """Returns the weblogo URI for its sites."""
        return bioutils.weblogo_uri(self.aligned_sites)


class MotifReport(MotifReportBase):
    """Collection of meta-sites that have the same TF and species."""
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


class EnsembleMotifReport(MotifReportBase):
    """Collection of meta-sites that are requested to be viewed together."""
    pass


def group_curation_site_instances(curation_site_instances):
    """Groups Curation_SiteInstance objects by TF-instance and genome."""
    TF_genome_pairs = curation_site_instances.filter(
        site_type='motif_associated').values_list(
        'site_instance__genome', 'curation__TF_instances').distinct()
    return [curation_site_instances.filter(site_instance__genome=genome,
                                           curation__TF_instances=TF_instances)
            for genome, TF_instances in TF_genome_pairs]


def build_motif_reports(curation_site_instances):
    """Given a list of Curation_SiteInstance objects, returns MotifReports."""
    curation_site_instance_grps = group_curation_site_instances(
        curation_site_instances)
    return [MotifReport(curation_site_instance_grp)
            for curation_site_instance_grp in curation_site_instance_grps]


def build_ensemble_report(curation_site_instances):
    """Given Curation_SiteInstance objects, returns the EnsembleMotifReport."""
    return EnsembleMotifReport(curation_site_instances)
