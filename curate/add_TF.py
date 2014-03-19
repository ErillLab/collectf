from forms import add_TF_form
from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages
from django.contrib.auth.decorators import user_passes_test
from django.contrib.auth.decorators import login_required
import models
import base

@login_required
def add_TF(request):
    """View for adding a TF"""
    if request.method == "POST":
        form = add_TF_form.AddTFForm(request.POST)
        if form.is_valid():
            cd = form.cleaned_data
            new_TF = models.TF(name=cd['name'],
                               description=cd['description'],
                               family=cd['family'])
            new_TF.save()
            messages.add_message(request, messages.INFO, "The TF was added successfully.")
            return HttpResponseRedirect(reverse(base.views.home))
    else:
        form = add_TF_form.AddTFForm()
        
    return render(request, "add_TF.html", {'form': form},
                  context_instance=RequestContext(request))

@login_required
def add_TF_family(request):
    """View for adding a TF family"""
    if request.method == "POST":
        form = add_TF_form.AddTFFamilyForm(request.POST)
        if form.is_valid():
            cd = form.cleaned_data
            new_family = models.TFFamily(name=cd['name'], description=cd['description'])
            new_family.save()
            messages.add_message(request, messages.INFO, "The TF family was added successfully.")
            return HttpResponseRedirect(reverse(base.views.home))
    else:
        form = add_TF_form.AddTFFamilyForm()

    return render(request, "add_TF_family.html", {'form': form},
                  context_instance=RequestContext(request))
