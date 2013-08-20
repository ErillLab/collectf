# Create your views here.
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect
from django.contrib import messages
from django.core.urlresolvers import reverse
import models
import random
import collectfapp.views
from baseapp import bioutils
from baseapp import lasagna
from django.core.mail import send_mail
import datetime
import logging

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

def feedback_send_email(request):
    your_email_address = request.POST['email']
    type = request.POST['type']
    comment = request.POST['comment']
    try:
        send_mail('CollecTF - %s feedback' % type,
                  comment,
                  your_email_address,
                  ['sefa1@umbc.edu', 'collectfdb@umbc.edu',],
                  fail_silently=False)
        messages.add_message(request, messages.INFO, 'Thanks for the feedback!')
        
    except:
        messages.add_message(request, messages.ERROR, "Something when wrong when sending feedback")

    return HttpResponseRedirect(reverse(collectfapp.views.home))

def register_request(request):
    if not request.POST:
        return render_to_response("register_form.html", {}, context_instance = RequestContext(request))
    # POST
    try:
        send_mail('CollecTF - new account request',
                  'username: %(username)s\nfirst_name: %(first_name)s\nlast_name: %(last_name)s\nemail: %(email)s' % request.POST,
                  request.POST['email'],
                  ['sefa1@umbc.edu', 'collectfdb@umbc.edu',],
                  fail_silently=False)
        messages.add_message(request, messages.INFO, 'Thanks for your request. We will contact you shortly for confirmation.')

    except:
        messages.add_message(request, messages.ERROR, "Something went wrong. Please try again.")

    return HttpResponseRedirect(reverse(collectfapp.views.home))
        

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


def log_IP(request):
    logger = logging.getLogger(__name__)
    logger.info('\t'.join(map(str, [datetime.datetime.now(),
                                    request.META['REMOTE_ADDR']))))
               
def greet(request):
    """Handler for main page right frame"""
    log_IP(request)
    template_file = "greet.html"
    random_rec = get_random_motif()
    return render_to_response(template_file,
                              {
                                  'random_rec': random_rec
                              },
                              context_instance=RequestContext(request))

def get_random_motif(motif_len_th=30, motif_sz_th=10):
    """Get random motif from the database, to display on the main page.
    Randomly selected motif must be no longer than <motif_len_th> and
    motif size must be at least <motif_sz_th>
    """

    TF_genome_list = models.Curation_SiteInstance.objects.values_list(
        'curation__TF__name',
        'curation__TF_instance__protein_accession',
        'site_instance__genome__genome_id',
        'site_instance__genome__genome_accession',
        'site_instance__genome__organism')

    while True:
        TF_genome = random.choice(TF_genome_list)
        # get all for this TF and genome
        csis = models.Curation_SiteInstance.objects.filter(curation__TF_instance__protein_accession=TF_genome[1],
                                                           site_instance__genome__genome_id=TF_genome[2],
                                                           is_motif_associated=True)
        csis = csis.all()
        if len(csis) > motif_sz_th:
            # align sites
            aligned, idxAligned, strands = lasagna.LASAGNA(map(lambda s:str(s.site_instance.seq).lower(), csis), 0)
            trimmed = lasagna.TrimAlignment(aligned) if len(aligned) > 1 else aligned
            trimmed = [s.upper() for s in trimmed]
            if len(trimmed[0]) < motif_len_th:
                # get non-motif-associated data
                ncsis = models.Curation_SiteInstance.objects.filter(curation__TF_instance__protein_accession=TF_genome[1],
                                                                   site_instance__genome__genome_id=TF_genome[2],
                                                                   is_motif_associated=False).all()

                return {
                    'aligned_sites': trimmed,
                    'TF_name': TF_genome[0],
                    'TF_accession': TF_genome[1],
                    'genome_accession': TF_genome[3],
                    'organism': TF_genome[4],
                    'view_all_csis': [csi.pk for csi in csis],
                    'view_all_ncsis': [ncsi.pk for ncsi in ncsis]
                    }
    
    
    
