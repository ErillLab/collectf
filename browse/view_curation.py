"""Handler to view curation"""
from base import models
from django.shortcuts import render_to_response
from django.template import RequestContext

def view_curation(request, cid):
    """Handler function for curation view"""
    curation = models.Curation.objects.get(curation_id=cid)
    return render_to_response("view_curation.html",
                              {'curation': curation},
                              context_instance = RequestContext(request))

