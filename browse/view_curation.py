"""Handler to view curation"""
from base import models
from django.shortcuts import render_to_response
from django.shortcuts import get_object_or_404
from django.template import RequestContext

def view_curation(request, cid):
    """Handler function for curation view"""
    curation = get_object_or_404(models.Curation, curation_id=cid)
    return render_to_response("view_curation.html",
                              {'curation': curation},
                              context_instance=RequestContext(request))

