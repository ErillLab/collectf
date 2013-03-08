from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
import models

@login_required
def view_curation(request, cid):
    """Handler function for curation view"""
    curation = models.Curation.objects.get(curation_id=cid)
    return HttpResponse("Hello view curation!")
