# Create your views here.
from django.shortcuts import render_to_response
from django.template import RequestContext

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

