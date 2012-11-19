"""
Publication submission views.
"""
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.contrib.formtools.preview import FormPreview
from pubform import PubmedPublicationForm
from pubform import NonPubmedPublicationForm
import bioutils
import sutils
import models
from models import Publication
import makeobj
import views
from django.contrib import messages

class PubSubmissionFormPreview(FormPreview):
    """Form preview view for publication submission"""
    FormPreview.form_template = "submit_pub.html"
    FormPreview.preview_template = "submit_pub_preview.html"
    
    # For process_preview and done methods, see Django docs
    def process_preview(self, request, form, context):
        cd = form.cleaned_data
        if 'pmid' in cd:
            # at this point, since pmid is validated in
            # pubform.PubmedPublicationForm.clean_pmid method, we're safe
            try:
                p = Publication.objects.get(pmid=cd['pmid'])
            except Publication.DoesNotExist:
                pubrec = bioutils.get_pubmed(cd['pmid'])
                p = makeobj.make_pub(pubrec, cd)
        else:  # non-pubmed publication
            pubrec = dict(Title=cd['title'],
                          AuthorList=cd['authors'].split(','),
                          FullJournalName=cd['journal'],
                          PubDate=cd['publication_date'],
                          Volume=cd['volume'],
                          Issue=cd['issue'],
                          Pages=cd['pages'])
            p = makeobj.make_pub(pubrec, cd)

        # put p to context to show as html on form preview
        context["pub"] = p
        sutils.sput(request.session, "pub_sub_publication", p)

    def done(self, request, cleaned_data):
        msg = "The paper was successfully submitted to be curated later."
        
        p = sutils.sget(request.session, "pub_sub_publication")
        if "assignment" in request.POST: # if assigned to the curator
            curator = models.Curator.objects.get(user=request.user)
            p.assigned_to = curator

        if "contains_no_data" in request.POST: # if the paper has no TFBS data
            note = " \nPaper has no TF-binding site data."
            p.submission_notes += note
            p.curation_complete = True
            msg = """The paper was marked as complete, since it does not have data"""

        p.save()  # insert into database

        messages.success(request, msg)
        
        return HttpResponseRedirect(reverse(views.home))

# pubmed publication submission handler
@login_required
def pubmed_submission(request):
    view = PubSubmissionFormPreview(PubmedPublicationForm)
    return view(request)

# non-pubmed publication submission handler
@login_required
def non_pubmed_submission(request):
    view = PubSubmissionFormPreview(NonPubmedPublicationForm)
    return view(request)
