"""
This module contains the class implementation for motif reports.
"""

import random

from django.utils.functional import cached_property

from core import bioutils
from core import lasagna
from core import metasite
from core import models
from core import motif_structure

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
    def species(self):
        """Returns the species of the sites in this motif report."""
        return self.meta_sites[0].delegate.site_instance.genome.taxonomy

    @property
    def species_name(self):
        """Returns the name of the species."""
        return self.meta_sites[0].delegate.site_instance.genome.organism

    @property
    def short_species_name(self):
        """Returns the species name in short format.

        For example, returns 'Y.pestis' for the species Yersinia pestis
        KIM10+.
        """
        genus, species = self.species_name.split()[:2]
        return genus[0] + '.' + species

    @property
    def genome_accession(self):
        """Returns the genome accession of the sites in the motif report."""
        return self.meta_sites[0].delegate.genome_accession

    @property
    def regulation_stats(self):
        """Returns the percentages of each regulation type."""
        stats = {}
        for tf_function in models.Curation_SiteInstance.TF_FUNCTION:
            sites = [meta_site for meta_site in self.meta_sites
                     if any(csi.TF_function == tf_function[0]
                            for csi in meta_site.curation_site_instances)]
            stats[tf_function[1]] = 100 * len(sites) / len(self.meta_sites)
        return stats

    @property
    def TF_conformation_stats(self):
        """Returns the percentages of TF-conformation."""
        stats = {}
        for tf_type in models.Curation_SiteInstance.TF_TYPE:
            sites = [meta_site for meta_site in self.meta_sites
                     if any(csi.TF_type == tf_type[0]
                            for csi in meta_site.curation_site_instances)]
            stats[tf_type[1]] = 100 * len(sites) / len(self.meta_sites)
        return stats

    @property
    def site_stats(self):
        """Returns the percentages of site types."""
        stats = {}
        for site_type in models.Curation_SiteInstance.SITE_TYPE:
            sites = [meta_site for meta_site in self.meta_sites
                     if any(csi.site_type == site_type[0]
                            for csi in meta_site.curation_site_instances)]
            stats[site_type[1]] = len(sites)
        return stats

    @property
    def GC_content(self):
        gc = sum(site.count('G') + site.count('C')
                 for meta_site in self.meta_sites
                 for site in meta_site.delegate_sequence)
        total = sum(len(site)
                    for meta_site in self.meta_sites
                    for site in meta_site.delegate_sequence)
        return float(gc) / total * 100

    @property
    def structure(self):
        """Returns the motif structure."""
        return motif_structure.find_pattern(self.aligned_sites)
    

class EnsembleMotifReport(MotifReportBase):
    """Collection of meta-sites that are requested to be viewed together."""
    pass


def group_curation_site_instances(curation_site_instances):
    """Groups Curation_SiteInstance objects by TF-instance and genome."""
    TF_genome_pairs = curation_site_instances.filter(
        site_type='motif_associated').order_by(
            'curation__TF_instances__TF',
            'site_instance__genome__organism').values_list(
                'site_instance__genome__organism',
                'curation__TF_instances').distinct()

    return [curation_site_instances.filter(
        site_instance__genome__organism=genome,
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


def random_motif_report(max_len=30, min_size=10):
    """Returns a random motif report.

    Randomly selects a motif until finds a good one.
    """
    TF_genome_pairs = models.Curation_SiteInstance.objects.values_list(
        'curation__TF_instances', 'site_instance__genome').distinct()
    while True:
        TF_instances, genome = random.choice(TF_genome_pairs)
        curation_site_instances = models.Curation_SiteInstance.objects.filter(
            curation__TF_instances=TF_instances,
            site_instance__genome=genome,
            site_type='motif_associated')
        if curation_site_instances.count() < min_size:
            continue
        reports = build_motif_reports(curation_site_instances)
        if len(reports[0].aligned_sites[0]) > max_len:
            continue
        return reports[0]
