import models

class MetaSite:
    """The idea behind meta-site is that a site may have been identified in
    multiple papers. To remove redundancy (i.e. displaying the same site
    sequence multiple times, one for each paper), they are collapsed into one
    "meta-site", which covers experimental evidence from all curations that the
    site is To be collapsed into the same meta-site, two site are required to
    overlap enough (but not necessarily 100%)."""

    def __init__(self, curation_site_instance):
        """Create a meta-site containing only one curation_site_instance."""
        self.cur_site_insts = []
        self.add_cur_site_inst(curation_site_instance)

    def membership_test(self, cur_site_inst):
        """Check if a curation-site-instance object is appropriate for the
        meta-site. Based on the type of cur_site_inst object, perform
        motif_associated_overlap_test or non_motif_associated_overlap_test"""
        if cur_site_inst.site_type == 'motif_associated':
            return (self.motif_associated_overlap_test(cur_site_inst) and
                    self.genome_test(cur_site_inst) and
                    self.TF_instances_test(cur_site_inst))
        elif cur_site_inst.site_type == 'non_motif_associated':
            return (self.non_motif_associated_overlap_test(cur_site_inst) and
                    self.genome_test(cur_site_inst) and
                    self.TF_instances_test(cur_site_inst))
        else:
            return False

    def TF_instances_test(self, cur_site_inst):
        """Check if two curation site instances have the same set of TF
        instances. Perform this after location overlap test since this takes
        more time and we want to eliminate as many candidates as possible via
        overlap test"""
        return set(self.TF_instances) == set(cur_site_inst.TF_instances)

    def genome_test(self, cur_site_inst):
        """Check if two curation site instances have the same genome"""
        return True # guaranteed to have same genome
        return self.genome == cur_site_inst.genome

    def motif_associated_overlap_test(self, cur_site_inst):
        """Given a curation-site-instance object, check if the meta-site and the
        curation-site-instance overlaps enough. The new site-instance is
        included into this meta-site if the overlap between new site and the
        delegate site is more than 75%"""
        def get_overlap(loca, locb):
            """Given two locations, return the overlap ratio."""
            overlap_len = max(0, min(loca[1], locb[1]) - max(loca[0], locb[0]))
            return float(overlap_len) / (loca[1]-loca[0]+1)


        loca = (self.delegate.site_instance.start, self.delegate.site_instance.end)
        locb = (cur_site_inst.site_instance.start, cur_site_inst.site_instance.end)
        overlap_a = get_overlap(loca, locb)
        overlap_b = get_overlap(locb, loca)
        return (overlap_a + overlap_b) / 2.0 >= 0.75

    def non_motif_associated_overlap_test(self, cur_site_inst):
        """Given a non-motif-associated curation-site-instance object, determine
        if it overlaps enough with this meta-site to be considered as a member
        of it. In case of enough overlap, the evidence from non-motif-associated
        site is integrated into the meta-site. The criteria is full overlap
        between the delegate-site (motif-associated one) and the target
        site-instance (non-motif_associated one)"""

        loca = (cur_site_inst.site_instance.start, cur_site_inst.site_instance.end)
        locb = (self.delegate.site_instance.start, self.delegate.site_instance.end)
        return (min(loca[0], loca[1]) <= min(locb[0], locb[1]) and
                max(loca[0], loca[1]) >= max(locb[0], locb[1]))

    def add_cur_site_inst(self, cur_site_inst):
        """Add a new curation-site-instance object to the meta-site"""
        self.cur_site_insts.append(cur_site_inst)
        if hasattr(self, '_techniques'):
            del self._techniques
        if hasattr(self, '_regulations'):
            del self._regulations
        if hasattr(self, '_curations'):
            del self._curations

    @property
    def genome(self):
        """Return genome of curation-site-instances."""
        if not hasattr(self, '_genome'):
            self._genome = self.cur_site_insts[0].genome
        return self._genome

    @property
    def TF_instances(self):
        """Return TF instances of curation-site-instances. After the
        curation-submission update on January 2014, a curation-site-instance can
        be linked to multiple TF-instance objects, to cover TFs composed of
        several subunits with different accession numbers (e.g. heterodimer such
        as IHF, composed of IHF alpha and beta)."""
        if not hasattr(self, '_TF_instances'):
           self._TF_instances =  self.cur_site_insts[0].TF_instances
        return self._TF_instances

    @property
    def delegate_sequence(self):
        if not hasattr(self, '_delegate_sequence'):
           self._delegate_sequence = self.cur_site_insts[0].site_instance.seq
        return self._delegate_sequence

    @property
    def delegate(self):
        """Each meta-site consists of a collection of curation-site-instances
        and the first one is used as a delegate."""
        return self.cur_site_insts[0]

    @property
    def delegate_site_instance(self):
        return self.delegate.site_instance

    @property
    def curations(self):
        """Return the list of curations for all curation-site-instance objects
        in the meta-site."""
        if not hasattr(self, '_curations'):
            self._curations = list(set([csi.curation for csi in self.cur_site_insts]))
        return self._curations

    @property
    def techniques(self):
        """Return the list of techniques that have been used to identify sites
        in the meta-site."""
        if not hasattr(self, '_techniques'):
            all_techniques = models.ExperimentalTechnique.objects.all()
            self._techniques = all_techniques.filter(curation_siteinstance__in=self.cur_site_insts)
        return self._techniques

    @property
    def regulations(self):
        """Collapse regulation information of all curation-site-instances into a
        single unit. For a single gene, there may be more than one regulation
        information."""
        if not hasattr(self, '_regulations'):
            regulations = models.Regulation.objects.filter(curation_site_instance__in=self.cur_site_insts)
            regs_exp_verified  = regulations.filter(evidence_type="exp_verified")
            regs_inferred      = regulations.filter(evidence_type="inferred")
            genes_exp_verified = regs_exp_verified.values_list("gene", flat=True)
            genes_inferred     = regs_inferred.values_list("gene", flat=True)
            self._regulations = ([reg for i,(gene,reg) in enumerate(zip(genes_exp_verified, regs_exp_verified))
                                  if gene not in genes_exp_verified[:i]] +
                                  [reg for i,(gene,reg) in enumerate(zip(genes_inferred, regs_inferred))
                                  if gene not in genes_inferred[:i] and gene not in genes_exp_verified])
        return self._regulations


def create_meta_sites(motif_cur_site_insts, non_motif_cur_site_insts):
    """Given a collection of motif-associated curation-site-instances and
    non-motif associated curation-site-instances, create meta-sites from
    motif-associated ones and integrate non-motif ones into them. In addition,
    integrate the experimental evidence, regulation information and curation
    information into the meta-site"""

    meta_sites = []
    # create meta-sites using motif-associated curation-site-instances
    for cur_site_inst in motif_cur_site_insts:
        # check if any meta-site overlaps enough
        for meta_site in meta_sites:
            if meta_site.membership_test(cur_site_inst):
                meta_site.add_cur_site_inst(cur_site_inst)
                break
        else: # It means none of the existing meta-sites are appropriate. Create
              # a new meta-site
            meta_sites.append(MetaSite(cur_site_inst))

    # integrate non-motif-associated curation-site-instances
    for cur_site_inst in non_motif_cur_site_insts:
        for meta_site in meta_sites:
            if meta_site.membership_test(cur_site_inst):
                meta_site.add_cur_site_inst(cur_site_inst)
                break
        else: # It means none of the existing meta-sites are appropriate. In the
              # case of non-motif-associated sites, DO NOT do anything.
              pass

    return meta_sites
