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


def citation(pubrec):
    """Create citation string from publication record"""
    title = pubrec["Title"]
    authors = pubrec["AuthorList"]
    journal = pubrec["FullJournalName"]
    return '|'.join([title, ','.join(authors), journal])

def make_pub(pubrec, cd):
    """Given publication record, retrieved from NCBI db (if pubmed publication),
    and entered user data, create models.Publication object and return it.
    """
    pmid = cd.get("pmid", None)  # None if not pubmed publication
    url ="http://www.ncbi.nlm.nih.gov/pubmed?term=%s" % pmid if pmid else cd['URL']
    publication_type = "pubmed" if pmid else "nonpubmed"
    p = Publication(publication_type = publication_type,
                    pmid=pmid,
                    citation=citation(pubrec),
                    url=url,
                    pdf=None,
                    contains_promoter_data=cd["contains_promoter_data"],
                    contains_expression_data=cd["contains_expression_data"],
                    submission_notes=cd["submission_notes"],
                    curation_complete=False)
    return p

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
                p = make_pub(pubrec, cd)
        else:  # non-pubmed publication
            pubrec = dict(Title=cd['title'],
                          AuthorList=cd['authors'].split(','),
                          FullJournalName=cd['journal'])
            p = make_pub(pubrec, cd)

        # put p to context to show as html on form preview
        context["pub"] = p


    def done(self, request, cleaned_data):
        p = request.session["publication"]
        p.save() # add pub to db
        if "assignment" in request.POST: # if assigned to the curator
            curator = models.Curator.objects.get(user=request.user)
            curator.assigned_papers.add(p)
        return HttpResponseRedirect("/success")

# pubmed publication submission handler
pubmed_submission = PubSubmissionFormPreview(PubmedPublicationForm)

# non-pubmed publication submission handler
non_pubmed_submission = PubSubmissionFormPreview(NonPubmedPublicationForm)

