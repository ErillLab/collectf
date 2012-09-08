from django import forms
import bioutils
import sutils
from models import Publication
                    
class PubmedPublicationForm(forms.Form):
    pmid = forms.CharField()
    contains_promoter_data = forms.BooleanField(required=False)
    contains_expression_data = forms.BooleanField(required=False)
    submission_notes = forms.CharField(widget=forms.Textarea, required=False)

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
    title = forms.CharField()
    authors = forms.CharField()
    journal = forms.CharField()
    URL = forms.URLField(verify_exists=False, required=False)
    # pdf?
    contains_promoter_data = forms.BooleanField(required=False)
    contains_expression_data = forms.BooleanField(required=False)
    submission_notes = forms.CharField(widget=forms.Textarea, required=False)
            
