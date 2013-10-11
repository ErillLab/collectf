# Create your views here.

import os
import pickle
from django.shortcuts import render
from django.template import RequestContext
from collectf import settings
import models

def curation_stats(request):
    """Handler for curation statistics page. Count the number of curations/sites for
    each TF and species in the database, pass dictionary to the HTML template"""
    pickle_file = os.path.join(settings.PICKLE_ROOT, 'dbstats.pickle')
    response_dict = pickle.load(open(pickle_file, 'r'))
    return render(request,
                  "database_stats.html",
                  response_dict,
                  context_instance=RequestContext(request))

def curator_roster(request):
    return render(request,
                  "curator_roster.html",
                  context_instance=RequestContext(request))


def experimental_techniques(request):
    exp_techniques = models.ExperimentalTechnique.objects.order_by('name').all()
    return render(request,
                  "experimental_techniques.html",
                  {'techs': exp_techniques},
                  context_instance=RequestContext(request))

    
    
def release_history(request):
    return render(request,
                  "release_history.html")
