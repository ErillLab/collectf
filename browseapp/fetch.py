"""Fetch objects from database"""

import models

def get_all_curations():
    return models.Curation.objects.all()


def get_all_TFs():
    return models.TF.objects.all()

def get_all_species():
    return models.Strain.objects.all()

def get_all_exp_techniques():
    return models.ExperimentalTechnique.objects.all()


