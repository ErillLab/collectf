"""
Publication submission views.
"""
from django.http import HttpResponseRedirect
from django.contrib.formtools.preview import FormPreview
from pubform import PubmedPublicationForm
import sutils
import models

class PubmedSubmissionFormPreview(FormPreview):
    """Form preview view for publication submission"""
    FormPreview.form_template = "submit_pub.html"
    FormPreview.preview_template = "submit_pub_preview.html"

    # workaround to access request.session
    def __call__(self, request, *a, **kw):
        self.form.request = request
        return super(PubmedSubmissionFormPreview, self).__call__(request, *a, **kw)
    
    # FOR process_preview and done methods, see Django docs
    def process_preview(self, request, form, context):
        context["pub"] = sutils.sget(request.session, "publication")

    def done(self, request, cleaned_data):
        p = request.session["publication"]
        p.save() # add pub to db
        if "assignment" in request.POST: # if assigned to the curator
            curator = models.Curator.objects.get(user=request.user)
            curator.assigned_papers.add(p)
        return HttpResponseRedirect("/success")

# pubmed publication submission handler
pubmed_submission = PubmedSubmissionFormPreview(PubmedPublicationForm)

