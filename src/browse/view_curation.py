"""Handler to view curation"""

from core import models
from core import metasite
from django.shortcuts import render
from django.shortcuts import get_object_or_404


def view_curation(request, curation_id):
    """Handler function for curation view."""
    curation = get_object_or_404(models.Curation, curation_id=curation_id)
    meta_sites = metasite.create_meta_sites(
        curation.curation_siteinstance_set.all())
    return render(request, 'view_curation.html', {'curation': curation,
                                                  'meta_sites': meta_sites})
