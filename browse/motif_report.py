"""
This file contains the class implementation for reports. A MotifReport object is
constructed from a collection of curation-site-instance objects that have the
same TF and species. For some uses of MotifReport class, all
curation-site-instance objects have to have exactly the same TF-instances
objects, but in other cases, having the same TF (there can be paralog TFs in the
same organism).
"""
import base
from base import metasite

class MotifReport:
    def __init__(self, m_cur_site_insts, nm_cur_site_insts=[]):
        """Given a collection of motif-associated and non-motif-associated
        curation site instances, create the report object. The constructor also
        checks if all curation_site_instance objects have the same TF and
        species"""
        # make sure that the list is not empty
        assert m_cur_site_insts
        self.m_cur_site_insts = m_cur_site_insts
        self.nm_cur_site_insts = nm_cur_site_insts

    def set_non_motif_curation_site_instances(self,non_motif_curation_site_insts):
        """Add some non-motif-associated curation-site-instances into the motif
        report object."""
        self.nm_cur_site_insts = non_motif_curation_site_insts

    def TF_check(self):
        """Check if all curation_site_instance objects have the same TF"""
        head = self.m_cur_site_insts[0]
        assert all(head.curation.TF == csi.curation.TF
                   for csi in self.m_cur_site_insts)
        assert all(head.curation.TF == csi.curation.TF
                   for csi in self.nm_cur_site_insts)

    def species_check(self):
        """Check if all curation_site_instance objects have the same species"""
        head = self.m_cur_site_insts[0]
        assert all(head.site_instance.genome.taxonomy == csi.site_instance.genome.taxonomy
                   for csi in self.m_cur_site_insts)
        assert all(head.site_instance.genome.taxonomy == csi.site_instance.genome.taxonomy
                   for csi in self.nm_cur_site_insts)

    def TF_accession_check(self):
        """Check if all curation-site-instances have the same TF-instance.
        Although all MotifReport objects have to satisfy the same TF and same
        species criteria, all TFs don't have to be exactly the same protein, or
        all species don't have to be same genomes. In other words, most of the
        cases, it is fine to run this function and get a False result."""
        head = self.m_cur_site_insts[0]
        return (all(head.curation.TF_instances == csi.curation.TF_instances
                    for csi in self.m_cur_site_insts) and
                all(head.curation.TF_instances == csi.curation.TF_instances
                    for csi in self.nm_cur_site_insts))
    
    def genome_accession_check(self):
        """Check if all curation-site-instances have the same
        genome-accession. Like TF_accession_check function, all
        curation-site-instance objects do not have to have the same genome."""
        head = self.m_cur_site_insts[0]
        return (all(head.site_instance.genome == csi.site_instance.genome
                    for csi in self.m_cur_site_insts) and
                all(head.site_instance.genome == csi.site_instance.genome
                    for csi in self.nm_cur_site_insts))

    @property
    def TF_name(self):
        """Return the name of the TF"""
        self.TF_accession_check
        return self.m_cur_site_insts[0].curation.TF.name

    @property
    def TF_accession(self):
        """Return the TF accession number.
        This function should be called only if all curation-site-instance
        objects have the same TF accession numbers."""
        get_TF = lambda csi: csi.curation.TF_instances.all()[0]
        assert all(get_TF(self.m_cur_site_insts[0]) == get_TF(x)
                   for x in self.m_cur_site_insts)
        return str(get_TF(self.m_cur_site_insts[0]).protein_accession)
    
    @property
    def species_name(self):
        """Return the name of the species"""
        return self.m_cur_site_insts[0].site_instance.genome.taxonomy.name

    @property
    def genome_accession(self):
        """Return the genome accession number.
        This function should be called only when all curation-site-instance
        objects have the same genome accession number."""
        get_genome = lambda csi: csi.site_instance.genome
        assert self.genome_accession_check()
        return str(get_genome(self.m_cur_site_insts[0]).genome_accession)
        
    def align_sites(self):
        """Align all binding sites using Lasagna"""
        return base.bioutils.run_lasagna([x.site_instance
                                          for x in self.m_cur_site_insts])

    def create_meta_sites(self):
        """Create meta-sites from curation-site-instances."""
        self.meta_sites = metasite.create_meta_sites(self.m_cur_site_insts,
                                                     self.nm_cur_site_insts)
        return self.meta_sites

    def get_all_motif_cur_site_insts(self,ids_only=True):
        """Return all motif-associated curation-site-instances (ids only by
        default; otherwise, return all objects)."""
        if ids_only:
            return [x.pk for x in self.m_cur_site_insts]
        else:
            return self.m_cur_site_insts

    @property
    def num_motif_cur_site_insts(self):
        """Return the number of motif-associated curation-site-instances."""
        return len(self.m_cur_site_insts)

    def get_all_non_motif_cur_site_insts(self,ids_only=True):
        """Return all non-motif-associated curation-site-instances (ids only by
        default; otherwise, return all objects)."""
        if ids_only:
            return [x.pk for x in self.nm_cur_site_insts]
        else:
            return self.nm_cur_site_insts

    def get_all_cur_site_insts(self, ids_only=True):
        """Return all curation-site-instance objects (both motif associated and
        non-motif associated)."""
        return (self.get_all_motif_cur_site_insts(ids_only) +
                self.get_all_non_motif_cur_site_insts(ids_only))

    def generate_browse_result_dict(self):
        """Generate a dictionary of values to add to the template context which
        will be rendered in browse-result page. The view will use these values
        just before rendering the template."""
        return {
            'TF_name': self.TF_name,
            'species_name': self.species_name,
            'cur_site_insts': self.get_all_cur_site_insts(),
        }

    def generate_view_reports_dict(self):
        """Generate a dictionary of values to add to the template context which
        will be rendered in view_reports page."""
        return {
            'TF_name': self.TF_name,
            'species_name': self.species_name,
            'meta_sites': self.create_meta_sites(),
            'aligned_sites': self.align_sites(),
            'cur_site_insts': self.get_all_cur_site_insts(),
        }
    
def make_reports(cur_site_insts):
    """Given a collection of motif-associated and non-motif-associated
    curation-site-instance objects, group them by TF and species and create
    MotifReport objects. If there is no motif-associated curation-site-instances for
    a given TF and species, a report object is not created, even if there are
    some non-motif-associated curation-site-instances. In other words, report
    objects are created for (TF,species) that have at least one motif-associated
    curation-site-instance."""
    # find all tuples of (TF,species)
    tf_species = cur_site_insts.values_list("curation__TF", "site_instance__genome__taxonomy")\
      .distinct()\
      .order_by("curation__TF__name", "site_instance__genome__taxonomy__name")
      
    # group all curation-site-instances by TF and species
    reports = []
    for (TF,species) in tf_species:
        # filter motif-assoc cur-site-insts and non-motif-assoc ones by TF and species
        m_csis = cur_site_insts.filter(curation__TF=TF,
                                       site_instance__genome__taxonomy=species,
                                       site_type="motif_associated")
        nm_csis = cur_site_insts.filter(curation__TF=TF,
                                        site_instance__genome__taxonomy=species,
                                        site_type="non_motif_associated")

        # generate a report only if there is at least one
        if m_csis: reports.append(MotifReport(m_csis, nm_csis))
    return reports

def make_ensemble_report(cur_site_insts):
    """Given a collection of curation-site-instance objects, align all of them
    and generate a report."""
    motif_associated = cur_site_insts.filter(site_type="motif_associated")
    non_motif_associated = cur_site_insts.filter(site_type="non_motif_associated")
    return MotifReport(motif_associated, non_motif_associated)
    

        


        
