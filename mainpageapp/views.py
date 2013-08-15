# Create your views here.
from django.shortcuts import render_to_response
from django.template import RequestContext
import models
import random

def about(request):
    template_file = "about.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def browse(request):
    template_file = "browse.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def search(request):
    template_file = "main_search.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def contribute(request):
    template_file = "contribute.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def feedback(request):
    template_file = "feedback.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def stats(request):
    template_file = "stats.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def cite(request):
    template_file = "cite.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def links(request):
    template_file = "links.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))

def acknowledgements(request):
    template_file = "acknowledgements.html"
    return render_to_response(template_file, {}, context_instance = RequestContext(request))


def greet(request):
    """Handler for main page right frame"""
    template_file = "sample.html"
    random_TF_instance, random_genome = get_random_motif()
    return render_to_response(template_file,
                              {
                                  'random_TF_instance': random_TF_instance,
                                  'random_genome': random_genome,
                              },
                              context_instance=RequestContext(request))

def get_random_motif(motif_len_th=30, motif_sz_th=10):
    """Get random motif from the database, to display on the main page.
    Randomly selected motif must be no longer than <motif_len_th> and
    motif size must be at least <motif_sz_th>
    """

    #select a random curation_site_instance
    random_id = random.randint(0, models.Curation_SiteInstance.objects.count()-1)
    random_csi = models.Curation_SiteInstance.objects.all()[random_id]

    return 1,2
