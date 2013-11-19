import add_TF_form
import views
import models
from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages

def add_TF(request):
    if request.method == "POST":
        form = add_TF_form.AddTFForm(request.POST)
        if form.is_valid():
            cd = form.cleaned_data
            new_TF = models.TF(name=cd['name'], description=cd['description'], family=cd['family'])
            messages.add_message(request, messages.INFO, "The TF was added.")
            return HttpResponseRedirect(reverse(views.home))

    else:
        form = add_TF_form.AddTFForm()
        
    return render(request, "add_TF.html",
                  {'form': form}, context_instance=RequestContext(request))
    
def add_TF_family(request):
    if request.method == "POST":
        form = add_TF_form.AddTFFamilyForm(request.POST)
        if form.is_valid():
            cd = form.cleaned_data
            new_family = models.TFFamily(name=cd['name'], description=cd['description'])
            new_family.save()
            messages.add_message(request, messages.INFO,
                                 "The TF family was added.")
            return HttpResponseRedirect(reverse(views.home))
    else:
        form = add_TF_form.AddTFFamilyForm()

    return render(request, "add_TF_family.html",
                  {'form': form}, context_instance=RequestContext(request))
