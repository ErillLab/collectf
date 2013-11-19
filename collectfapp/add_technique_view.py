import add_technique_form
from django.shortcuts import render
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.core.urlresolvers import reverse
from django.contrib import messages
import models
import views

def add_technique(request):
    if request.method == "POST":
        form = add_technique_form.AddTechniqueForm(request.POST)
        if form.is_valid():
            cd = form.cleaned_data
            new_technique = models.ExperimentalTechnique(name=cd['name'],
                                                         description=cd['description'],
                                                         preset_function=cd['type'])
            new_technique.save()
            for cat in cd['categories']:
                new_technique.categories.add(cat)
            messages.add_message(request, messages.INFO, "The experimental technique was added.")
    else:
        form = add_technique_form.AddTechniqueForm()

    return render(request, "add_technique.html",
                  {'form': form}, context_instance=RequestContext(request))
