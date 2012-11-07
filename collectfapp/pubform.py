from django import forms
import bioutils
import sutils
from models import Publication
                    
class PubmedPublicationForm(forms.Form):
    pmid = forms.CharField(label="PMID", required=True)
    reported_TF = forms.CharField(label="Reported TF(s)", required=False)
    reported_species = forms.CharField(label="Reported species", required=False)

    contains_promoter_data = forms.BooleanField(required=False,
                                                label="The manuscript contains promoter information")
    contains_expression_data = forms.BooleanField(required=False,
                                                  label="The manuscript contains expression information")
    submission_notes = forms.CharField(widget=forms.Textarea, required=False,
                                       label="Submission notes")

    def clean_pmid(self):
        cp = self.cleaned_data['pmid'] # cleaned pmid
        # check if in database
        try:
            Publication.objects.get(pmid=cp)
            raise forms.ValidationError("Publication is already in database")
        except Publication.DoesNotExist:
            pass
        
        if not bioutils.get_pubmed(cp):
            raise forms.ValidationError("Invalid PMID")
        return cp

class NonPubmedPublicationForm(forms.Form):
    # form for non-pubmed publications or unpublished data
    authors = forms.CharField(label="Authors (Name Initials,)")
    title = forms.CharField(label="Title")
    journal = forms.CharField(label="Journal")
    year = forms.CharField(label="Year")
    volume = forms.CharField(label="Volume")
    pages = forms.CharField(label="Pages (xxx-yyy)")
    URL = forms.URLField(verify_exists=False, required=False, label="URL")
    # pdf?
    reported_TF = forms.CharField(label="Reported TF(s)", required=False)
    reported_species = forms.CharField(label="Reported species",
                                       required=False)

    contains_promoter_data = forms.BooleanField(required=False,
                                                label="The manuscript contains promoter information")
    contains_expression_data = forms.BooleanField(required=False,
                                                  label="The manuscript contains expression information")
    submission_notes = forms.CharField(widget=forms.Textarea, required=False,
                                       label="Submission notes")
            
