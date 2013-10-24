from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect

def home(request):
    template_file = "main.html"
    return render_to_response(template_file, {}, context_instance=RequestContext(request))
