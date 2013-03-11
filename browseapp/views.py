from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.template import RequestContext
import fetch

def view_all_curations(request):
    """Handler function to see all curations at once"""
    template_vals = {
        "curations": fetch.get_all_curations()
    }
    return render_to_response("curation_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))


def browse_get(request):
    template_vals = {
        'TFs': fetch.get_all_TFs(),
        'species': fetch.get_all_species(),
        'experimental_techniques': fetch.get_all_exp_techniques(),
    }
    return render_to_response("browse.html",
                              template_vals,
                              context_instance=RequestContext(request))

def browse_post(request):
    pass

def browse(request):
    """Handler function to browse database"""
    if not request.POST:
        return browse_get(request)

    return browse_post(request)
    
    

