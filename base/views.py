"""Main views"""

from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.core.urlresolvers import reverse
import django.contrib.auth.views
import homepage.views


def home(request):
    if request.user.is_authenticated():
        return home_user(request)
    else:
        return home_guest(request)

def home_user(request):
    return HttpResponseRedirect(reverse(homepage.views.greet))

def home_guest(request):
    return HttpResponseRedirect(reverse(homepage.views.greet))

# registration handler
def register(request, redirect=''):
    """View function for registering user. The user is linked with
    models.Curator object."""
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

def get_weblogo(request):
    """Given a POST request having a list of binding sites, return the generated
    weblogo."""
    sites = request.POST['sites'].strip().split(',')
    sites = map(str, sites) # make them string
    # make sure they all have the same length
    assert all(len(sites[0]) == len(site) for site in sites)
    return HttpResponse(bioutils.weblogo_uri(sites), "text/plain")
