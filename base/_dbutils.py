"""This file contains implementations for some frequently performed database
interactions."""

def get_all_cur_site_insts():
    """Return all curation-site-instance objects."""
    return models.Curation_SiteInstance.objects.all()

def get_all_motif_assoc_cur_site_insts():
    """Return all motif associated curation-site-inst objects"""
    return get_all_cur_site_insts().filter(site_type="motif_associated")

def get_all_non_motif_assoc_cur_site_insts():
    """Return all non-motif associated curation-site-inst objects"""
    return get_all_cur_site_insts().filter(site_type="non_motif_associated")
