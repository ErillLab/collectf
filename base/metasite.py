from base import models

class MetaSite:
    """MetaSite class definition.

    A site could be reported in multiple papers. To avoid redundancy
    (i.e. presenting the same site sequence multiple times, one per paper), the
    Curation_SiteInstance objects are collapsed into one entity called a
    "meta-site". To be collapsed into the same meta-site, two sites are required
    to overlap more than a defined threshold.
    """

    def __init__(self, curation_site_instance):
        self.curation_site_instances = [curation_site_instance]

    def add(self, curation_site_instance):
        """Adds the given Curation_SiteInstance object to the meta-site."""
        self.curation_site_instances.append(curation_site_instance)

    @property
    def delegate(self):
        """Returns the delegate SiteInstance object."""
        return self.curation_site_instances[0]
        
    @property
    def genome_accession(self):
        """Returns the genome accession of the meta-site."""
        return self.delegate.site_instance.genome.genome_accession

    @property
    def TF_instance_accessions(self):
        """Returns the accession numbers of TF_instances of the meta-site."""
        return self.delegate.curation.TF_instances.values_list(
            'uniprot_accession', flat=True)

    @property
    def motif_id(self):
        """Returns the motif ID of its Curation_SiteInstances.

        A TF may have multiple motifs for the same species and this field is
        used to identify the specific motif that the binding site belongs to.
        """
        return self.delegate.motif_id

    @property
    def site_type(self):
        """Returns the site type of the meta-site.

        The site type is the site type of the delegate site.
        """
        return self.delegate.site_type

    @property
    def techniques(self):
        """Returns the experimental techniques for all Curation_SiteInstances"""
        technique_ids = set(
            technique.technique_id
            for curation_site_instance in self.curation_site_instances
            for technique in curation_site_instance.experimental_techniques.all())
        print technique_ids
        return models.ExperimentalTechnique.objects.filter(
            technique_id__in=technique_ids)

    @property
    def regulations(self):
        regulations = models.Regulation.objects.filter(
            curation_site_instance__in=self.curation_site_instances)
        exp_verified_regulations = regulations.filter(
            evidence_type='exp_verified')
        exp_verified_genes = exp_verified_regulations.values_list(
            'gene', flat=True)
        inferred_regulations = regulations.filter(evidence_type='inferred')
        inferred_genes = inferred_regulations.values_list('gene', flat=True)
        meta_site_regulations = (
            [regulation
             for i, (gene, regulation) in
             enumerate(zip(exp_verified_genes, exp_verified_regulations))
             if gene not in exp_verified_genes[:i]] +
            [regulation
             for i, (gene, regulation)
             in enumerate(zip(inferred_genes, inferred_regulations))
             if gene not in exp_verified_genes and
             gene not in inferred_genes[:i]])
        return meta_site_regulations

    @property
    def curation_ids(self):
        """Returns the curation IDs of the meta-site."""
        return [curation_site_instance.curation.curation_id
                for curation_site_instance in self.curation_site_instances]

    @property
    def curation_site_instance_ids(self):
        """Returns the Curation_SiteInstance IDs of the meta-site."""
        return [curation_site_instance.pk
                for curation_site_instance in self.curation_site_instances]

    def membership_test(self, curation_site_instance):
        """Checks if curation_site_instance can be member of the meta-site.

        Based on the type of Curation_SiteInstance object, performs
        motif_associated_overlap_test or non_motif_associated_overlap_test.
        """
        if curation_site_instance.site_type in ['motif_associated',
                                                'var_motif_associated']:
            return (
                self.genome_test(curation_site_instance) and
                self.TF_instances_test(curation_site_instance) and
                self.motif_id_test(curation_site_instance) and
                self.motif_associated_overlap_test(curation_site_instance))
        elif cur_site_inst.site_type == 'non_motif_associated':
            return (
                self.genome_test(curation_site_instance) and
                self.TF_instances_test(curation_site_instance) and
                self.non_motif_associated_overlap_test(
                    curation_site_instance))

        return False

    def TF_instances_test(self, curation_site_instance):
        """Checks if a curation_site_instance have the same TF instances."""
        return (self.delegate.curation.TF_instances ==
                curation_site_instance.curation.TF_instances)

    def genome_test(self, curation_site_instance):
        """Checks if the curation_site_instance the same genome."""
        return (self.genome_accession ==
                curation_site_instance.site_instance.genome.genome_accession)

    def motif_id_test(self, curation_site_instance):
        """Checks if the curation_site_instance is from the same motif."""
        return self.motif_id == curation_site_instance.motif_id

    def motif_associated_overlap_test(self, curation_site_instance):
        """Checks if the meta-site and the Curation_SiteInstance overlaps.

        Includes the new Curation_SiteInstance into this meta-site if the
        overlap between new site and the delegate site is more than 75%.
        """
        def get_overlap(loca, locb):
            """Given two locations, returns the overlap ratio."""
            overlap_len = max(0, min(loca[1], locb[1]) - max(loca[0], locb[0]))
            return float(overlap_len) / (loca[1]-loca[0]+1)

        loca = (self.delegate.site_instance.start,
                self.delegate.site_instance.end)
        locb = (curation_site_instance.site_instance.start,
                curation_site_instance.site_instance.end)
        overlap_a = get_overlap(loca, locb)
        overlap_b = get_overlap(locb, loca)
        return (overlap_a + overlap_b) / 2.0 >= 0.75

    def non_motif_associated_overlap_test(self, curation_site_instance):
        """Checks if the meta-site and the Curation_SiteInstance overlaps.

        In case of enough overlap, integrates the evidence from
        non-motif-associated site into the meta-site. The criteria is full
        overlap between the delegate-site (motif-associated one) and the target
        site-instance (non-motif_associated one)
        """
        loca = (curation_site_instance.site_instance.start,
                curation_site_instance.site_instance.end)
        locb = (self.delegate.site_instance.start,
                self.delegate.site_instance.end)
        return (min(loca[0], loca[1]) <= min(locb[0], locb[1]) and
                max(loca[0], loca[1]) >= max(locb[0], locb[1]))

def create_meta_sites(curation_site_instances):
    """Creates meta sites from given Curation_SiteInstance objects.

    Given a collection of motif-associated Curation_SiteInstances and non-motif
    associated Curation_SiteInstances, creates meta-sites from motif-associated
    ones and integrates non-motif ones into them.
    """
    meta_sites = []
    # Create meta-sites using motif-associated Curation_SiteInstances.
    motif_associated = curation_site_instances.filter(
        site_type='motif_associated')
    for curation_site_instance in motif_associated:
        # Check if any meta-site overlaps enough
        for meta_site in meta_sites:
            if meta_site.membership_test(curation_site_instance):
                meta_site.add(curation_site_instance)
                break
        else:
            # None of the existing meta-sites are good. Create a new meta-site.
            meta_sites.append(MetaSite(curation_site_instance))

    variable_motif_associated = curation_site_instances.filter(
        site_type='var_motif_associated')
    for curation_site_instance in variable_motif_associated:
        for meta_site in meta_sites:
            if meta_site.membership_test(curation_site_instance):
                meta_site.add(curation_site_instance)
                break
        else:
            meta_sites.append(MetaSite(curation_site_instance))

    # Integrate non-motif-associated curation-site-instances
    non_motif_associated = curation_site_instances.filter(
        site_type='non_motif_associated')
    for curation_site_instance in non_motif_associated:
        for meta_site in meta_sites:
            if meta_site.membership_test(curation_site_instance):
                meta_site.add(curation_site_instance)
                break
        else:
            # None of the existing meta-sites are good. In the case of
            # non-motif-associated sites, DO NOT do anything.
            pass

    return meta_sites
