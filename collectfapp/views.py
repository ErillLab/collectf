from django.shortcuts import render_to_response
from django.contrib.auth.decorators import login_required
import models
from django.template import RequestContext


# Create your views here.

# home view
#@login_required
def home(request):
    """Home view function to choose an action.
    The user can submit paper for curation, curate a paper or edit one of their
    own previous curations.
    """
    template_vals = dict()
    if request.user.is_authenticated():
        curator,created = models.Curator.objects.get_or_create(user=request.user)
        curations = curator.curation_set.all()
        template_vals["user"] = request.user
        template_vals["curator"] = curator
        template_vals["curations"] = curations
        
    return render_to_response("choose.html", template_vals,
                              context_instance = RequestContext(request))

@login_required
def success(request):
    """Success view handler"""
    template_vals = {"user": request.user}
    return render_to_response("success.html", template_vals,
                              context_instance = RequestContext(request))
