from django.shortcuts import render_to_response
from django.contrib.auth.decorators import login_required
from django.http import HttpResponseRedirect
import models
from django.template import RequestContext
from django.core.mail import send_mail
from django.contrib import messages
import mainpageapp.views
from django.core.urlresolvers import reverse

# Create your views here.

# home view
#@login_required
def home(request):
    """Home view function to choose an action.
    The user can submit paper for curation, curate a paper or edit one of their
    own previous curations.
    """
    template_vals = dict()
    if request.user.is_authenticated():
        curator,created = models.Curator.objects.get_or_create(user=request.user)
        curations = curator.curation_set.all()
        template_vals["user"] = request.user
        template_vals["curator"] = curator
        template_vals["curations"] = curations
        template_file = "choose.html"
        return render_to_response(template_file, template_vals,
                                  context_instance = RequestContext(request))
    
    else:
        #template_file = "greet.html"
        print 'goto greet'
        return HttpResponseRedirect(reverse(mainpageapp.views.greet))
        


@login_required
def success(request):
    """Success view handler"""
    template_vals = {"user": request.user}
    return render_to_response("success.html", template_vals,
                              context_instance = RequestContext(request))



def pub_external_submission(request):
    if not request.POST:
        return render_to_response("pub_external_submission.html", {},
                                  context_instance = RequestContext(request))
    # POST
    try:
        send_mail('CollecTF - new paper submission',
                  'email: %s, pmid: %s' % (request.POST['email'], request.POST['pmid']),
                  request.POST['email'],
                  ['sefa1@umbc.edu', 'collectfdb@umbc.edu',],
                  fail_silently=False
                  )
        messages.add_message(request, messages.INFO, 'Thanks for your submission. We will include it in our database soon.')

    except:
        messages.add_message(request, messages.ERROR, "Something went wrong. Please try again.")

    return home(request)
