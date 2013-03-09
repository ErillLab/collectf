"""Fetch objects from database"""

import models

def get_all_curations():
    return models.Curation.objects.all()


