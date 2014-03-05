"""View functions for publication submission"""
from django.http import HttpResponseRedirect
from django.contrib.auth.decorators import login_required
from django.core.urlresolvers import reverse
from django.contrib.formtools.preview import FormPreview
from forms.add_publication_form import PubmedPublicationForm
from forms.add_publication_form import NonPubmedPublicationForm
from django.contrib import messages
import session_utils
import base
import create_object
from models import Publication
from models import Curator

class PubSubmissionFormPreview(FormPreview):
    """Form preview view for publication submission. """
    FormPreview.form_template = "add_publication.html"
    FormPreview.preview_template = "add_publication_preview.html"
    
    # For process_preview and done methods, see Django docs
    def process_preview(self, request, form, context):
        cd = form.cleaned_data
        if 'pmid' in cd:
            # At this point, since pmid is validated in
            # pubform.PubmedPublicationForm.clean_pmid method, we're safe
            try:
                # check if publication is already in the database.
                p = Publication.objects.get(pmid=cd['pmid'])
            except Publication.DoesNotExist:
                # if not, add to the DB.
                pubrec = base.bioutils.get_pubmed(cd['pmid'])
                p = create_object.make_pub(pubrec, cd)
        else:  # non-pubmed publication, create it manually
            pubrec = dict(Title=cd['title'],
                          AuthorList=cd['authors'].split(','),
                          FullJournalName=cd['journal'],
                          PubDate=cd['publication_date'],
                          Volume=cd['volume'],
                          Issue=cd['issue'],
                          Pages=cd['pages'])
            p = create_object.make_pub(pubrec, cd)

        # At this point, the Publication object is created, BUT it is not in DB yet.
        # Pass the object for review.
        context["pub"] = p
        session_utils.put(request.session, "publication", p)

    def done(self, request, cleaned_data):
        """Add object to the database and return the appropriate message to the user."""
        msg = "The paper was successfully submitted to be curated later."
        p = session_utils.get(request.session, "publication")
        if "assignment" in request.POST: # if assigned to the curator
            curator,_ = Curator.objects.get_or_create(user=request.user)
            p.assigned_to = curator

        if "contains_no_data" in request.POST: # if the paper has no TFBS data
            note = " \nPaper has no TF-binding site data."
            p.submission_notes += note
            p.curation_complete = True
            msg = """The paper was marked as complete, since it does not have data."""

        p.save()  # insert into database
        messages.success(request, msg)
        return HttpResponseRedirect(reverse(base.views.home))


# pubmed publication submission handler
@login_required
def pubmed_submission(request):
    """Handler for paper submission that has a PMID."""
    view = PubSubmissionFormPreview(PubmedPublicationForm)
    return view(request)

# non-pubmed publication submission handler
@login_required
def non_pubmed_submission(request):
    """Handler for paper submission that has not a PMID."""
    view = PubSubmissionFormPreview(NonPubmedPublicationForm)
    return view(request)
