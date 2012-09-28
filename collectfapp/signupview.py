import django.contrib.auth.views
from django.contrib import auth
from django.http import HttpResponseRedirect
from django.shortcuts import render_to_response
from signupform import *
from django.template import RequestContext
from models import Curator

# registration handler
def register(request, redirect=''):
    """View function for registering user. The user is linked with
    models.Curator object.
    """
    if request.method == "GET":
        form = CuratorRegistrationForm()
        
    if request.method == "POST":
        form = CuratorRegistrationForm(request.POST)
        if form.is_valid():
            new_user = form.save()           
            Curator(user=new_user).save() # create Curator object in db
            new_user = auth.authenticate(username=request.POST.get("username"),
                                         password=request.POST.get("password1"))
            auth.login(request, new_user)
            return HttpResponseRedirect(redirect)

    return render_to_response("registration/register.html", {"form": form},
                              context_instance=RequestContext(request))

# let django handle login
login = django.contrib.auth.views.login
# and logout
logout =  django.contrib.auth.views.logout

