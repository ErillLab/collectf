from django.contrib import messages
from django.shortcuts import get_object_or_404
from django.shortcuts import render

from core import dbxref
from core import models
from core import metasite


def view_site(request, dbxref_id):
    """Handler to view site instances."""
    curation_site_instance = get_object_or_404(
        models.Curation_SiteInstance, pk=dbxref.from_ncbi_dbxref(dbxref_id))
    if curation_site_instance.is_obsolete:
        messages.add_message(
            request, messages.ERROR,
            "The requested site instance seems obsolete. "
            "It will be removed on the next release. "
            "Description: %s" % curation_site_instance.why_obsolete)

    meta_site = metasite.MetaSite(curation_site_instance)

    return render(request, 'view_site.html',
                  {'dbxref': dbxref_id,
                   'meta_site': meta_site})
