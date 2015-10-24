import models

class MetaSite:
    """Definition for CollecTF meta-sites.
    
    The idea behind meta-site is that a site may have been identified in
    multiple papers. To remove redundancy (i.e. displaying the same site
    sequence multiple times, one for each paper), they are collapsed into one
    "meta-site", which covers experimental evidence from all curations that the
    site is To be collapsed into the same meta-site, two site are required to
    overlap enough (but not necessarily 100%).
    """

    def __init__(self, curation_site_instance):
        """Creates a meta-site containing only one Curation_SiteInstance."""
        self.cur_site_insts = []
        self.add_cur_site_inst(curation_site_instance)

    def membership_test(self, cur_site_inst):
        """Checks if a Curation_SiteInstance object is good for the meta-site.
        
        Based on the type of Curation_SiteInstance object, performs
        motif_associated_overlap_test or non_motif_associated_overlap_test.
        """
        if cur_site_inst.site_type in ['motif_associated',
                                       'var_motif_associated']:
            return (self.motif_associated_overlap_test(cur_site_inst) and
                    self.genome_test(cur_site_inst) and
                    self.TF_instances_test(cur_site_inst) and
                    self.motif_id_test(cur_site_inst))
        elif cur_site_inst.site_type == 'non_motif_associated':
            return (self.non_motif_associated_overlap_test(cur_site_inst) and
                    self.genome_test(cur_site_inst) and
                    self.TF_instances_test(cur_site_inst))
        else:
            return False

    def TF_instances_test(self, cur_site_inst):
        """Checks if two Curation_SiteInstances have the same TF instances.

        Performs this after location overlap test since this takes more time and
        we want to eliminate as many candidates as possible via overlap test.
        """
        return set(self.TF_instances) == set(cur_site_inst.TF_instances)

    def genome_test(self, cur_site_inst):
        """Checks if two Curation_SiteInstances have the same genome."""
        return (self.cur_site_insts[0].site_instance.genome ==
                cur_site_inst.site_instance.genome)
    
    def motif_id_test(self, cur_site_inst):
        """Checks if two site instances belong to the same motif."""
        return self.motif_id == cur_site_inst.motif_id

    def motif_associated_overlap_test(self, cur_site_inst):
        """Checks if the meta-site and the Curation_SiteInstance overlaps.

        Includes the new Curation_SiteInstance into this meta-site if the
        overlap between new site and the delegate site is more than 75%.
        """
        def get_overlap(loca, locb):
            """Given two locations, returns the overlap ratio."""
            overlap_len = max(0, min(loca[1], locb[1]) - max(loca[0], locb[0]))
            return float(overlap_len) / (loca[1]-loca[0]+1)

        loca = (self.delegate.site_instance.start, self.delegate.site_instance.end)
        locb = (cur_site_inst.site_instance.start, cur_site_inst.site_instance.end)
        overlap_a = get_overlap(loca, locb)
        overlap_b = get_overlap(locb, loca)
        return (overlap_a + overlap_b) / 2.0 >= 0.75

    def non_motif_associated_overlap_test(self, cur_site_inst):
        """Checks if the meta-site and the Curation_SiteInstance overlaps.

        In case of enough overlap, integrates the evidence from
        non-motif-associated site into the meta-site. The criteria is full
        overlap between the delegate-site (motif-associated one) and the target
        site-instance (non-motif_associated one)
        """
        loca = (cur_site_inst.site_instance.start,
                cur_site_inst.site_instance.end)
        locb = (self.delegate.site_instance.start,
                self.delegate.site_instance.end)
        return (min(loca[0], loca[1]) <= min(locb[0], locb[1]) and
                max(loca[0], loca[1]) >= max(locb[0], locb[1]))

    def add_cur_site_inst(self, cur_site_inst):
        """Adds a new Curation_SiteInstance object to the meta-site."""
        self.cur_site_insts.append(cur_site_inst)

    @property
    def genome_accession(self):
        """Returns the genome accession of the meta-site."""
        return self.cur_site_insts[0].site_instance.genome.genome_accession

    @property
    def TF_instances(self):
        """Returns TF instances of the Curation_SiteInstances."""
        return self.cur_site_insts[0].TF_instances

    @property
    def motif_id(self):
        """Returns the motif ID of its Curation_SiteInstances.

        A TF may have multiple motifs for the same species and this field is
        used to identify the specific motif that the binding site belongs to.
        """
        return self.cur_site_insts[0].motif_id

    @property
    def site_type(self):
        """Returns the site type of the meta-site.

        The site type is the site type of the delegate site.
        """
        return self.cur_site_insts[0].site_type

    @property
    def delegate_sequence(self):
        """Returns the delegate sequence of the meta-site."""
        return self.cur_site_insts[0].site_instance.seq

    @property
    def delegate(self):
        """Returns the delegate Curation_SiteInstance object.

        Each meta-site consists of a collection of curation-site-instances and
        the first one is used as a delegate.
        """
        return self.cur_site_insts[0]

    @property
    def delegate_site_instance(self):
        """Returns the delegate SiteInstance object."""
        return self.delegate.site_instance

    @property
    def curations(self):
        """Returns the list of curations for all Curation_SiteInstances."""
        return list(set(csi.curation for csi in self.cur_site_insts))

    @property
    def techniques(self):
        """Returns the list of techniques that have been used."""
        all_techniques = models.ExperimentalTechnique.objects.all()
        return all_techniques.filter(
            curation_siteinstance__in=self.cur_site_insts)

    @property
    def regulations(self):
        """Merges regulation information of all Curation_SiteInstances."""
        regulations = models.Regulation.objects.filter(
            curation_site_instance__in=self.cur_site_insts)
        regs_exp_verified = regulations.filter(evidence_type='exp_verified')
        regs_inferred = regulations.filter(evidence_type='inferred')
        genes_exp_verified = regs_exp_verified.values_list('gene', flat=True)
        genes_inferred = regs_inferred.values_list("gene", flat=True)
        # TODO(sefa): simplify this.
        return ([reg for i, (gene,reg) in enumerate(zip(genes_exp_verified, regs_exp_verified))
                 if gene not in genes_exp_verified[:i]] +
                [reg for i,(gene,reg) in enumerate(zip(genes_inferred, regs_inferred))
                 if gene not in genes_inferred[:i] and gene not in genes_exp_verified])

def create_meta_sites(motif_cur_site_insts,
                      non_motif_cur_site_insts,
                      var_motif_cur_site_insts=None):
    """Creates meta sites from given Curation_SiteInstance objects.

    Given a collection of motif-associated Curation_SiteInstances and non-motif
    associated Curation_SiteInstances, creates meta-sites from motif-associated
    ones and integrates non-motif ones into them. In addition, integrate the
    experimental evidence, regulation information and curation information into
    the meta-site.
    """

    meta_sites = []
    # Create meta-sites using motif-associated Curation_SiteInstances.
    for cur_site_inst in motif_cur_site_insts:
        # Check if any meta-site overlaps enough
        for meta_site in meta_sites:
            if meta_site.membership_test(cur_site_inst):
                meta_site.add_cur_site_inst(cur_site_inst)
                break
        else:
            # None of the existing meta-sites are good. Create a new meta-site.
            meta_sites.append(MetaSite(cur_site_inst))

    if var_motif_cur_site_insts:
        for cur_site_inst in var_motif_cur_site_insts:
            for meta_site in meta_sites:
                if meta_site.membership_test(cur_site_inst):
                    meta_site.add_cur_site_inst(cur_site_inst)
                    break
            else:
                meta_sites.append(MetaSite(cur_site_inst))

    # Integrate non-motif-associated curation-site-instances
    for cur_site_inst in non_motif_cur_site_insts:
        for meta_site in meta_sites:
            if meta_site.membership_test(cur_site_inst):
                meta_site.add_cur_site_inst(cur_site_inst)
                break
        else:
            # None of the existing meta-sites are good. In the case of
            # non-motif-associated sites, DO NOT do anything.
            pass

    return meta_sites
