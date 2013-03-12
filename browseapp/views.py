from django.http import HttpResponse
from django.shortcuts import render_to_response
from django.shortcuts import render
from django.template import RequestContext

import fetch
import forms

def view_all_curations(request):
    """Handler function to see all curations at once"""
    template_vals = {
        "curations": fetch.get_all_curations()
    }
    return render_to_response("curation_view_all.html",
                              template_vals,
                              context_instance=RequestContext(request))


def browse_get(request):
    form = forms.BrowseForm()
    return render(request,
                  "browse.html",
                  {'form': forms.BrowseForm(),},
                  context_instance=RequestContext(request))

def browse_post(request):
    TF = request.POST.get('TF')
    species = request.POST.get('species')
    csi = fetch.get_curations(TF, species)
    response_dict = {
        'form': forms.BrowseForm(initial={'TF': TF, 'species': species}),
        'csi': csi if csi else False,
    }
    return render(request,
                  "browse.html",
                  response_dict,
                  context_instance=RequestContext(request))

    

def browse(request):
    """Handler function to browse database"""
    if not request.POST:
        return browse_get(request)

    return browse_post(request)
    
    

