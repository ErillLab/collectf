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
    TF_id = request.POST.get('TF')
    species_id = request.POST.get('species')
    csi = fetch.get_curations(TF_id, species_id)
    response_dict = {
        'form': forms.BrowseForm(initial={'TF': TF_id, 'species': species_id}),
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

def report_FASTA(request, TF_id, species_id):
    """Given a list of sites, report FASTA file containing sites for particular TF
    and species"""
    # get TF and species from id
    TF = fetch.get_TF_by_id(TF_id)
    species = fetch.get_species_by_id(species_id) 

    # set HttpResponse stuff
    response = HttpResponse(content_type='application/download')
    response['Content-Disposition'] = 'attachment;filename="%s_%s.fasta"' % (TF.name, species.name)
    
    # write all sites to file
    for curation_site_instance in fetch.get_curations(TF_id, species_id):
        response.write(curation_site_instance.site_instance.to_fasta())
        
    return response
    
