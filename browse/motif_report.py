"""
This file contains the class implementation for reports. A MotifReport object is
constructed from a collection of meta-sites that have the same TF and species.
"""

from base import models
from base import metasite
from base import bioutils

class MotifReport:
    def __init__(self, curation_site_instances):
        """Creates a MotifReport objcet.

        Given a collection of Curation_SiteInstance objects, creates the report
        object.
        """
        self.meta_sites = list(set(
            curation_site_instance.meta_site
            for curation_site_instance in curation_site_instances))

    @property
    def id(self):
        return id(self)

    @property
    def TF(self):
        """Returns the TF."""
        return self.meta_sites[0].delegate.curation.TF

    @property
    def TF_name(self):
        """Returns the name of the TF."""
        return self.TF.name

    @property
    def TF_accessions(self):
        """Returns the TF accession numbers."""
        return self.meta_sites[0].TF_instance_accessions

    @property
    def species_name(self):
        """Returns the name of the species."""
        return self.meta_sites[0].delegate.site_instance.genome.organism

    @property
    def species(self):
        """Returns the species."""
        return self.meta_sites[0].delegate.site_instance.genome.taxonomy

    @property
    def genome_accession(self):
        """Returns the genome accession number."""
        return self.meta_sites[0].genome_accession

    @property
    def curation_site_instance_ids(self):
        """Returns the Curation_SiteInstance object ids."""
        return [csi for meta_site in self.meta_sites
                for csi in meta_site.curation_site_instance_ids]

    def align_sites(self):
        """Aligns all binding sites using Lasagna."""
        site_instances = [meta_site.delegate.site_instance
                          for meta_site in self.meta_sites
                          if meta_site.site_type == 'motif_associated']
        return bioutils.run_lasagna(site_instances)

def make_reports(curation_site_instances):
    """Generates motif reports from given Curation_SiteInstance objects.

    Given a collection of motif-associated and non-motif-associated
    Curation_SiteInstance objects, groups them by TF and species and creates
    MotifReport objects.

    If there is no motif-associated Curation_SiteInstance for a given TF and
    species, doesn't create a report object, even if there are some
    non-motif-associated Curation_SiteInstance objects.
    """
    # Find all tuples of (TF, species)
    tf_species = curation_site_instances.values_list(
        'curation__TF_instances__TF',
        'site_instance__genome__taxonomy').distinct().order_by(
            'curation__TF_instances__TF__name',
            'site_instance__genome__taxonomy__name')

    # Group all curation-site-instances by TF and species
    reports = []
    for (TF, species) in tf_species:
        filtered_group = curation_site_instances.filter(
            curation__TF_instances__TF=TF,
            site_instance__genome__taxonomy=species,
            site_type='motif_associated')
        if filtered_group:
            reports.append(MotifReport(filtered_group))

    return reports

def make_ensemble_report(curation_site_instances):
    """Aligns all Curation_SiteInstance objects and generates a report."""
    motif_associated_curation_site_instances = curation_site_instances.filter(
        site_type='motif_associated')
    if motif_associated_curation_site_instances:
        return MotifReport(motif_associated_curation_site_instances)

# TODO(sefa): New class for ensemble reports
