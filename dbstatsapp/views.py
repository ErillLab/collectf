# Create your views here.

import pickle
from django.shortcuts import render
from django.template import RequestContext

def curation_stats(request):
    """Handler for curation statistics page. Count the number of curations/sites for
    each TF and species in the database, pass dictionary to the HTML template"""
    pickle_file = 'dbstats.pickle'
    response_dict = pickle.load(open(pickle_file, 'r'))
    return render(request,
                  "database_stats.html",
                  response_dict,
                  context_instance=RequestContext(request))
