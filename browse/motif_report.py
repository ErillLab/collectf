"""
This file contains the class implementation for reports. A MotifReport object is
constructed from a collection of curation-site-instance objects that have the
same TF and species. For some uses of MotifReport class, all
curation-site-instance objects have to have exactly the same TF-instances
objects, but in other cases, having the same TF (there can be paralog TFs in the
same organism).
"""

from tqdm import tqdm

import base
from base import bioutils
from base import metasite
import base.models as models

class MotifReport:
    def __init__(
            self, m_cur_site_insts,
            nm_cur_site_insts=models.Curation_SiteInstance.objects.none(),
            var_cur_site_insts= models.Curation_SiteInstance.objects.none()):
        """Creates a MotifReport objcet.
        
        Given a collection of motif-associated and non-motif-associated
        Curation_SiteInstance objects, creates the report object. Also checks if
        all Curation_SiteInstance objects have the same TF and species.
        """
        # Make sure that the list is not empty
        self.m_cur_site_insts = m_cur_site_insts
        self.nm_cur_site_insts = nm_cur_site_insts
        self.var_cur_site_insts = var_cur_site_insts

        # Compute and store additional data
        self.compute_meta_sites()
        self.align_meta_sites()
        self.build_sequence_logo()

    def TF_check(self):
        """Checks if all Curation_SiteInstance objects have the same TF."""
        head = self.m_cur_site_insts[0]
        assert all(head.curation.TF == csi.curation.TF
                   for csi in self.m_cur_site_insts)
        assert all(head.curation.TF == csi.curation.TF
                   for csi in self.nm_cur_site_insts)
        assert all(head.curation.TF == csi.curation.TF
                   for csi in self.var_cur_site_insts)

    def species_check(self):
        """Checks if all Curation_SiteInstance objects have the same species."""
        head = self.m_cur_site_insts[0]
        assert all(head.site_instance.genome.taxonomy ==
                   csi.site_instance.genome.taxonomy
                   for csi in self.m_cur_site_insts)
        assert all(head.site_instance.genome.taxonomy ==
                   csi.site_instance.genome.taxonomy
                   for csi in self.nm_cur_site_insts)
        assert all(head.site_instance.genome.taxonomy ==
                   csi.site_instance.genome.taxonomy
                   for csi in self.var_cur_site_insts)

    def TF_accession_check(self):
        """Checks if all Curation_SiteInstances  have the same TF-instance.
        
        Although all MotifReport objects have to satisfy the same TF and same
        species criteria, all TFs don't have to be exactly the same protein, or
        all species don't have to be same genomes. In other words, most of the
        cases, it is fine to run this function and get a False result.
        """
        head = self.m_cur_site_insts[0]
        return (all(head.curation.TF_instances == csi.curation.TF_instances
                    for csi in self.m_cur_site_insts) and
                all(head.curation.TF_instances == csi.curation.TF_instances
                    for csi in self.nm_cur_site_insts) and
                all(head.curation.TF_instances == csi.curation.TF_instances
                    for csi in self.var_cur_site_insts))

    def genome_accession_check(self):
        """Checks if all Curation_SiteInstances have the same genome accession.

        Like TF_accession_check function, all Curation_SiteInstance objects do
        not have to have the same genome.
        """
        head = self.m_cur_site_insts[0]
        return (all(head.site_instance.genome == csi.site_instance.genome
                    for csi in self.m_cur_site_insts) and
                all(head.site_instance.genome == csi.site_instance.genome
                    for csi in self.nm_cur_site_insts) and
                all(head.site_instance.genome == csi.site_instance.genome
                    for csi in self.var_cur_site_insts))

    @property
    def TF_name(self):
        """Returns the name of the TF."""
        #self.TF_accession_check()
        return self.m_cur_site_insts[0].curation.TF

    @property
    def TF(self):
        """Returns the TF."""
        return self.m_cur_site_insts[0].curation.TF

    @property
    def TF_accessions(self):
        """Returns the TF accession number.
        
        Should be called only if all Curation_SiteInstance objects have the
        same TF accession numbers.
        """
        get_TF = lambda csi: csi.curation.TF_instances.all()
        return [str(TF_inst.uniprot_accession)
                for TF_inst in get_TF(self.m_cur_site_insts[0])]

    @property
    def species_name(self):
        """Returns the name of the species."""
        return self.m_cur_site_insts[0].site_instance.genome.organism

    @property
    def species(self):
        """Returns the species."""
        return self.m_cur_site_insts[0].site_instance.genome.taxonomy

    @property
    def genome_accession(self):
        """Returns the genome accession number.
        
        Should be called only when all Curation_SiteInstance objects have the
        same genome accession number.
        """
        get_genome = lambda csi: csi.site_instance.genome
        return str(get_genome(self.m_cur_site_insts[0]).genome_accession)

    def compute_meta_sites(self):
        """Creates meta-sites from curation-site-instances."""
        self.meta_sites = metasite.create_meta_sites(self.m_cur_site_insts,
                                                     self.nm_cur_site_insts,
                                                     self.var_cur_site_insts)
    def get_meta_sites(self):
        return self.meta_sites

    def align_meta_sites(self):
        """Aligns all binding sites using Lasagna."""
        meta_sites = self.get_meta_sites()
        self.aligned_meta_sites = bioutils.run_lasagna(
            [x.delegate_site_instance for x in meta_sites
             if x.site_type == 'motif_associated'])

    def get_aligned_meta_sites(self):
        return self.aligned_meta_sites

    def build_sequence_logo(self):
        seqs = self.get_aligned_meta_sites()
        self.weblogo_uri = bioutils.weblogo_uri(seqs)

    @property
    def num_motif_cur_site_insts(self):
        """Returns the number of motif-associated Curation_SiteInstances."""
        return len(self.m_cur_site_insts)

    def get_all_motif_cur_site_insts(self):
        """Returns all motif-associated CurationSiteInstances."""
        return self.m_cur_site_insts.all()

    def get_all_motif_cur_site_insts_ids(self):
        """Returns IDs of all motif-associated Curation_SiteInstances."""
        return [x.pk for x in self.m_cur_site_insts]

    def get_all_non_motif_cur_site_insts(self):
        """Returns all non-motif-associated Curation_SiteInstances."""
        return self.nm_cur_site_insts.all()

    def get_all_non_motif_cur_site_insts_ids(self):
        """Returns IDs of all non-motif-associated Curation_SiteInstances."""
        return [x.pk for x in self.nm_cur_site_insts]

    def get_all_var_motif_cur_site_insts(self):
        """Returns all variable-motif-associated Curation_SiteInstances."""
        return self.var_cur_site_insts.all()

    def get_all_motif_cur_site_insts_ids(self):
        """Return IDs of variable-motif-associated Curation_SiteInstances."""
        return [x.pk for x in self.var_cur_site_insts]

    def get_all_cur_site_insts(self):
        """Returns all Curation_SiteInstance object."""
        return (self.get_all_motif_cur_site_insts() |
                self.get_all_non_motif_cur_site_insts() |
                self.get_all_var_motif_cur_site_insts())

    def get_all_cur_site_insts_ids(self):
        """Returns IDs of all Curation_SiteInstance objects."""
        return [x.pk for x in self.get_all_cur_site_insts()]

    def generate_browse_result_dict(self):
        """Generates a dictionary of values to render browse result template."""
        return {
            'TF': self.TF,
            'species': self.species,
        }

    def generate_view_reports_dict(self):
        """Generates a dictionary of values to render view result template."""
        meta_site_info = [
            {'delegate_sequence': ms.delegate_sequence,
             'delegate_site_instance': ms.delegate_site_instance,
             'genome_accession': ms.genome_accession,
             'site_type': ms.site_type,
             'TF_instances': ms.TF_instances,
             'motif_id': ms.motif_id,
             'techniques': ms.techniques,
             'regulations': ms.regulations,
             'curations': ms.curations}
            for ms in self.get_meta_sites()]
        return {
            'TF_name': self.TF_name,
            'species_name': self.species_name,
            'TF_accessions': self.TF_accessions,
            'genome_accession': self.genome_accession,
            'id': id(self),
            'meta_sites': meta_site_info,
            'aligned_sites': self.get_aligned_meta_sites(),
            'weblogo': self.weblogo_uri,
            'cur_site_insts': self.get_all_cur_site_insts_ids(),
        }

def make_reports(cur_site_insts):
    """Generates motif reports from given Curation_SiteInstance objects.
    
    Given a collection of motif-associated and non-motif-associated
    Curation_SiteInstance objects, groups them by TF and species and creates
    MotifReport objects. If there is no motif-associated Curation_SiteInstance
    for a given TF and species, doesn't create a report object, even if there
    are some non-motif-associated Curation_SiteInstance objects.
    """
    # Find all tuples of (TF, species)
    tf_species = cur_site_insts.values_list(
        "curation__TF_instances__TF", "site_instance__genome__taxonomy")\
            .distinct()\
            .order_by("curation__TF_instances__TF__name",
                      "site_instance__genome__taxonomy__name")

    # Group all curation-site-instances by TF and species
    reports = []
    for (TF, species) in tf_species:
        m_csis = cur_site_insts.filter(
            curation__TF_instances__TF=TF,
            site_instance__genome__taxonomy=species,
            site_type='motif_associated')
        nm_csis = cur_site_insts.filter(
            curation__TF_instances__TF=TF,
            site_instance__genome__taxonomy=species,
            site_type='non_motif_associated')
        var_csis = cur_site_insts.filter(
            curation__TF_instances__TF=TF,
            site_instance__genome__taxonomy=species,
            site_type='var_motif_associated')

        # Generate a report only if there is at least one
        if m_csis:
            reports.append(MotifReport(m_csis, nm_csis, var_csis))

    return reports

def make_distinct_reports(cur_site_insts):
    """Generates motif reports by TF instances.

    Similar to make_reports. The difference is that rather than grouping
    curation_site_instances by TF and species names, groups them by accession
    numbers.
    """
    tf_species = cur_site_insts.values_list(
        'curation__TF_instances__TF_instance_id').distinct()

    reports = []
    for tf_insts in tf_species:
        m_csis = cur_site_insts.filter(curation__TF_instances=tf_insts,
                                       site_type='motif_associated')
        nm_csis = cur_site_insts.filter(curation__TF_instances=tf_insts,
                                        site_type='non_motif_associated')
        var_csis = cur_site_insts.filter(curation__TF_instances=tf_insts,
                                         site_type='var_motif_associated')
            

        # Generate only if there is at least one motif-associated site.
        if m_csis:
            reports.append(MotifReport(m_csis, nm_csis, var_csis))

    return reports

def make_ensemble_report(cur_site_insts):
    """Aligns all Curation_SiteInstance objects and generates a report."""
    motif_associated = cur_site_insts.filter(site_type='motif_associated')
    non_motif_associated = cur_site_insts.filter(
        site_type='non_motif_associated')
    var_motif_associated = cur_site_insts.filter(
        site_type='var_motif_associated')
    if motif_associated:
        return MotifReport(motif_associated, non_motif_associated,
                           var_motif_associated)
