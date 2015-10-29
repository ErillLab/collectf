from tqdm import tqdm

from base import models

def clear_meta_sites():
    """Clears meta-sites from the database."""
    for curation_site_instance in tqdm(models.Curation_SiteInstance.objects.all()):
        curation_site_instance.meta_site = None
        curation_site_instance.save()

    models.MetaSite.objects.all().delete()

def compute_meta_sites_for_motif_associated_sites(curation_site_instances):
    """Creates meta-sites for the given Curation_SiteInstance objects."""
    for curation_site_instance in tqdm(curation_site_instances):
        if curation_site_instance.meta_site:
            continue
        for meta_site in models.MetaSite.objects.all():
            if meta_site.membership_test(curation_site_instance):
                curation_site_instance.meta_site = meta_site
                curation_site_instance.save()
                break
        else:
            # Create a new meta-site.
            meta_site = models.MetaSite(delegate=curation_site_instance)
            meta_site.save()
            curation_site_instance.meta_site = meta_site
            curation_site_instance.save()

def compute_meta_sites_for_non_motif_associated_sites(curation_site_instances):
    """Creates meta-sites for the given non-motif-associated
    Curation_SiteInstance objects."""
    for curation_site_instance in tqdm(curation_site_instances):
        for meta_site in models.MetaSite.objects.all():
            if meta_site.membership_test(curation_site_instance):
                curation_site_instance.meta_site = meta_site
                curation_site_instance.save()
        else:
            # Do not create meta-site with leading non-motif-associated
            # Curation_SiteInstance objects.
            pass

def compute_meta_site_regulations():
    """Computes regulated genes for all meta-sites."""
    for meta_site in tqdm(models.MetaSite.objects.all()):
        regulations = models.Regulation.objects.filter(
            curation_site_instance__in=meta_site.curation_siteinstance_set.all())
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
        for regulation in meta_site_regulations:
            regulation.meta_site = meta_site
            regulation.save()

def run():
    # # Clear meta-sites
    # clear_meta_sites()
    
    # Motif-associated Curation_SiteInstances
    compute_meta_sites_for_motif_associated_sites(
        models.Curation_SiteInstance.objects.filter(
            site_type='motif_associated'))
    # Variable motif-associated Curation_SiteInstances
    compute_meta_sites_for_motif_associated_sites(
        models.Curation_SiteInstance.objects.filter(
            site_type='var_motif_associated'))
    # Non-motif-associated Curation_SiteInstances
    compute_meta_sites_for_non_motif_associated_sites(
        models.Curation_SiteInstance.objects.filter(
            site_type='none_motif_associated'))

    # Compute meta-site regulations
    compute_meta_site_regulations()
