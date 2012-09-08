"""
Publication submission views.
"""
from django.http import HttpResponseRedirect
from django.contrib.formtools.preview import FormPreview
from pubform import PubmedPublicationForm
from pubform import NonPubmedPublicationForm
import bioutils
import sutils
import models
from models import Publication
import makeobj

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
                          FullJournalName=cd['journal'])
            p = makeobj.make_pub(pubrec, cd)

        # put p to context to show as html on form preview
        context["pub"] = p
        sutils.sput(request.session, "publication", p)

    def done(self, request, cleaned_data):
        p = sutils.sget(request.session, "publication")
        p.save() # add pub to db
        if "assignment" in request.POST: # if assigned to the curator
            curator = models.Curator.objects.get(user=request.user)
            curator.assigned_papers.add(p)
        return HttpResponseRedirect("/success")

# pubmed publication submission handler
pubmed_submission = PubSubmissionFormPreview(PubmedPublicationForm)

# non-pubmed publication submission handler
non_pubmed_submission = PubSubmissionFormPreview(NonPubmedPublicationForm)

