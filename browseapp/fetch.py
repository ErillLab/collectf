"""Fetch objects from database"""

import models
import baseapp.bioutils as bioutils

def get_all_curations():
    return models.Curation.objects.all()

def get_curation_site_instances(TF, species):
    """Return all Curation_SiteInstance objects for particular TF and species"""
    curation_site_instances = models.Curation_SiteInstance.objects.filter(
        site_instance__genome__strain=species,
        curation__TF=TF)
    return curation_site_instances

def get_all_publications():
    return models.Publication.objects.all()

def get_all_TFs():
    return models.TF.objects.order_by('name')

def get_TFs_by_family(TF_family):
    return models.TF.objects.filter(family=TF_family).order_by('name')

def get_TF_family_by_id(TF_family_id):
    return models.TFFamily.objects.get(pk=TF_family_id)

def get_all_TF_families():
    return models.TFFamily.objects.order_by('name')

def get_all_species():
    return models.Strain.objects.order_by('name')

def get_lineage(species_taxon_id):
    return bioutils.get_taxon_info(species_taxon_id)

def get_all_exp_techniques():
    return models.ExperimentalTechnique.objects.all()
    
def get_all_TF_species():
    return models.TF.objects.values('strain')

def get_all_site_instances():
    return models.SiteInstance.objects.all()
    
def get_site_instances_by_species(strain):
    return models.SiteInstance.objects.filter(genome__strain=strain)

def get_TF_by_id(TF_id):
    return models.TF.objects.get(TF_id=TF_id)

def get_species_by_id(species_id):
    return models.Strain.objects.get(pk=species_id)

