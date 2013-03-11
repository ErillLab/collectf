from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
import models

@login_required
def view_curation(request, cid):
    """Handler function for curation view"""
    curation = models.Curation.objects.get(curation_id=cid)
    template_vals = dict(curation=curation)
    
    return render_to_response("curation_view.html", template_vals,
                              context_instance = RequestContext(request))


