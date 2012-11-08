from django import forms
import bioutils
import sutils
from models import Publication
                    
class PubmedPublicationForm(forms.Form):
    pmid = forms.CharField(label="PMID", required=True,
                           help_text="""Paste the Pubmed ID number obtained from
                           the NCBI website.""") 
    reported_TF = forms.CharField(label="Reported TF(s)", required=False,
                                  help_text="""Type the name of the
                                  transcription factor(s) reported in the
                                  manuscript.""")
    reported_species = forms.CharField(label="Reported species", required=False,
                                       help_text="""Type the name of the species
                                       reported in the manuscript.""")

    contains_promoter_data = forms.BooleanField(required=False,
            label="The manuscript contains promoter information",
            help_text="""The paper provides experimental data on the structure
            and sequence of TF-regulated promoter.""")
    contains_expression_data = forms.BooleanField(required=False,
                label="The manuscript contains expression information",
                help_text="""The paper provides experimental support for
            TF-mediated regulation of genes.""")
    submission_notes = forms.CharField(widget=forms.Textarea, required=False,
                                       label="Submission notes",
                                       help_text="""Include any additional
details about the submission. For instance, you might indicate the approximate
number of sites reported, whether high-throughput techniques are used or any
other factor that might help prioritize curation.""")

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
    publication_date = forms.CharField(label="Publication date")
    volume = forms.CharField(label="Volume")
    issue = forms.CharField(label="Issue")
    pages = forms.CharField(label="Pages (xxx-yyy)")
    URL = forms.URLField(verify_exists=False, required=False, label="URL")
    # pdf?
    
    reported_TF = forms.CharField(label="Reported TF(s)", required=False,
                                  help_text="""Type the name of the
                                  transcription factor(s) reported in the
                                  manuscript.""")
    reported_species = forms.CharField(label="Reported species", required=False,
                                       help_text="""Type the name of the species
                                       reported in the manuscript.""")

    contains_promoter_data = forms.BooleanField(required=False,
            label="The manuscript contains promoter information",
            help_text="""The paper provides experimental data on the structure
            and sequence of TF-regulated promoter.""")
    contains_expression_data = forms.BooleanField(required=False,
                label="The manuscript contains expression information",
                help_text="""The paper provides experimental support for
            TF-mediated regulation of genes.""")
    submission_notes = forms.CharField(widget=forms.Textarea, required=False,
                                       label="Submission notes",
                                       help_text="""Include any additional
details about the submission. For instance, you might indicate the approximate
number of sites reported, whether high-throughput techniques are used or any
other factor that might help prioritize curation.""")
            
