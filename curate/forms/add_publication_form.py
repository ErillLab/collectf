from django import forms
from base import bioutils
from curate import session_utils
from curate.models import Publication
import help_texts

ht = help_texts.pubmed_publication_form
class PubmedPublicationForm(forms.Form):
    pmid = forms.CharField(label="PMID", required=True, help_text=ht['pmid'])
    
    reported_TF = forms.CharField(label="Reported TF(s)", required=False, help_text=ht['reported_TF'])
    
    reported_species = forms.CharField(label="Reported species", required=False, help_text=ht['reported_species'])
    
    contains_promoter_data = forms.BooleanField(required=False,
            label="The manuscript contains promoter information",
            help_text=ht['contains_promoter_data'])
    
    contains_expression_data = forms.BooleanField(required=False,
                label="The manuscript contains expression information",
                help_text=ht['contains_expression_data'])
    
    submission_notes = forms.CharField(widget=forms.Textarea,
                                       required=False,
                                       label="Submission notes",
                                       help_text=ht['submission_notes'])

    def clean_pmid(self):
        cp = self.cleaned_data['pmid'].strip() # cleaned pmid
        # check if in database
        try:
            Publication.objects.get(pmid=cp)
            raise forms.ValidationError("Publication is already in database")
        except Publication.DoesNotExist:
            pass
        
        if not bioutils.get_pubmed(cp):
            raise forms.ValidationError("Invalid PMID")
        return cp


ht = help_texts.pubmed_publication_form
class NonPubmedPublicationForm(forms.Form):
    # form for non-pubmed publications or unpublished data
    authors = forms.CharField(label="Authors (Name Initials,)")
    title = forms.CharField(label="Title")
    journal = forms.CharField(label="Journal")
    publication_date = forms.CharField(label="Publication date")
    volume = forms.CharField(label="Volume")
    issue = forms.CharField(label="Issue")
    pages = forms.CharField(label="Pages (xxx-yyy)")
    URL = forms.URLField(required=False, label="URL")

    reported_TF = forms.CharField(label="Reported TF(s)", required=False,
                                  help_text=ht['reported_TF'])

    reported_species = forms.CharField(label="Reported species", required=False,
                                       help_text=ht['reported_species'])

    contains_promoter_data = forms.BooleanField(required=False,
            label="The manuscript contains promoter information",
            help_text=ht["contains_promoter_data"])
    
    contains_expression_data = forms.BooleanField(required=False,
                label="The manuscript contains expression information",
                help_text=ht["contains_expression_data"])
    
    submission_notes = forms.CharField(widget=forms.Textarea, required=False,
                                       label="Submission notes",
                                       help_text=ht["submission_notes"])
            
