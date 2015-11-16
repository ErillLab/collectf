"""New publication FormPreview definition using django-formtools."""

from django.contrib.auth.decorators import login_required
from django.contrib import messages
from django.shortcuts import redirect


from formtools.preview import FormPreview

from .forms.add_publication import PubMedPublicationForm
from .forms.add_publication import NonPubMedPublicationForm
from core import entrez_utils
from core import models


class AddPublicationFormPreview(FormPreview):
    """Form preview view for publication submission."""
    FormPreview.form_template = "add_publication.html"
    FormPreview.preview_template = "add_publication_preview.html"
    
    def process_preview(self, request, form, context):
        cleaned_data = form.cleaned_data
        if 'pmid' in cleaned_data:
            # PubMed publication
            record = entrez_utils.get_pubmed(cleaned_data['pmid'])
            publication = models.Publication(
                publication_type='pubmed',
                pmid=cleaned_data['pmid'],
                authors=', '.join(record.get('AuthorList')),
                journal=record.get('FullJournalName'),
                title=unicode(record.get('Title')),
                publication_date=record.get('PubDate'),
                volume=record.get('Volume'),
                issue=record.get('Issue'),
                pages=record.get('Pages'),
                url=cleaned_data.get('url'),
                pdf=None,
                contains_promoter_data=cleaned_data['contains_promoter_data'],
                contains_expression_data=cleaned_data['contains_expression_data'],
                submission_notes=cleaned_data['submission_notes'],
                curation_complete=False,
                reported_TF=cleaned_data['reported_TF'],
                reported_species=cleaned_data['reported_species'])
        else:
            # Non-pubmed publication, create it manually
            publication = models.Publication(
                publication_type='nonpubmed',
                authors=cleaned_data['authors'],
                journal=cleaned_data['journal'],
                title=cleaned_data['title'],
                publication_date=cleaned_data['publication_date'],
                volume=cleaned_data['volume'],
                issue=cleaned_data['issue'],
                pages=cleaned_data['pages'],
                url=cleaned_data.get('url'),
                pdf=None,
                contains_promoter_data=cleaned_data['contains_promoter_data'],
                contains_expression_data=cleaned_data['contains_expression_data'],
                submission_notes=cleaned_data['submission_notes'],
                curation_complete=False,
                reported_TF=cleaned_data['reported_TF'],
                reported_species=cleaned_data['reported_species'])
            
        # Store the Publication object for review.
        context['publication'] = publication
        request.session['publication'] = publication

    def done(self, request, cleaned_data):
        """Adds Publication object to the database."""    
        publication =  request.session['publication']
        # Check if the paper is assigned to the curator.
        if 'assignment' in request.POST: 
            curator, _ = models.Curator.objects.get_or_create(
                user=request.user)
            publication.assigned_to = curator

        # Check if the paper is marked as not having data.
        if 'contains_no_data' in request.POST:
            publication.submission_notes += "\nPaper has no TF-binding site data."
            publication.curation_complete = True
            msg = """The paper was marked as complete, since it does not have
            data."""
        else:
            msg = "The paper was successfully submitted to be curated later."

        publication.save()  # save to the database.
        messages.success(request, msg)
        return redirect('homepage_home')

    
@login_required
def add_pubmed_publication(request):
    """View function for adding a new PubMed publication."""
    return AddPublicationFormPreview(PubMedPublicationForm)(request)


@login_required
def add_non_pubmed_publication(request):
    """View function for adding a new non-PubMed publication."""
    return AddPublicationFormPreview(NonPubMedPublicationForm)(request)
