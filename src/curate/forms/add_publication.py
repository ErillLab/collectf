from django import forms

from core import models
from core import entrez_utils

class PubMedPublicationForm(forms.Form):
    pmid = forms.CharField(
        label="PMID",
        required=True,
        help_text="Enter the PubMed ID obtained from the NCBI website.")

    reported_TF = forms.CharField(
        label="Reported TF(s)",
        required=False,
        help_text="""Type the name of the transcription factor(s) reported in
        the manuscript.""")

    reported_species = forms.CharField(
        label="Reported species",
        required=False,
        help_text="Type the name of the species reported in the manuscript.")

    contains_promoter_data = forms.BooleanField(
        required=False,
        label="The manuscript contains promoter information",
        help_text="""The paper provides experimental data on the structure and
        sequence of TF-regulated promoter.""")

    contains_expression_data = forms.BooleanField(
        required=False,
        label="The manuscript contains expression information",
        help_text="""The paper provides experimental support for TF-mediated
        regulation of genes.""")
    
    submission_notes = forms.CharField(
        label="Submission notes",
        required=False, 
        help_text="""Include any additional details about the submission. For
        instance, you might indicate the approximate number of sites reported,
        whether high-throughput techniques are used or any other factor that
        might help prioritize curation.""",
        widget=forms.Textarea)


    def clean_pmid(self):
        """Cleans PubMed ID field."""
        pmid = self.cleaned_data['pmid'].strip()
        try:
            models.Publication.objects.get(pmid=pmid)
            raise forms.ValidationError("Publication is already in database.")
        except models.Publication.DoesNotExist:
            try:
                entrez_utils.get_pubmed(pmid)
            except entrez_utils.EntrezException:
                raise forms.ValidationError("Invalid PubMed ID.")
        return pmid


class NonPubMedPublicationForm(forms.Form):
    authors = forms.CharField(label="Authors (Name Initials,)")
    
    title = forms.CharField(label="Title")
    
    journal = forms.CharField(label="Journal")
    
    publication_date = forms.CharField(label="Publication date")
    
    volume = forms.CharField(label="Volume")
    
    issue = forms.CharField(label="Issue")
    
    pages = forms.CharField(label="Pages (xxx-yyy)")
    
    URL = forms.URLField(label="URL",
                         required=False)

    reported_TF = forms.CharField(
        label="Reported TF(s)",
        required=False,
        help_text="""Type the name of the transcription factor(s) reported in
        the manuscript.""")

    reported_species = forms.CharField(
        label="Reported species",
        required=False,
        help_text="Type the name of the species reported in the manuscript.")

    contains_promoter_data = forms.BooleanField(
        required=False,
        label="The manuscript contains promoter information",
        help_text="""The paper provides experimental data on the structure and
        sequence of TF-regulated promoter.""")

    contains_expression_data = forms.BooleanField(
        required=False,
        label="The manuscript contains expression information",
        help_text="""The paper provides experimental support for TF-mediated
        regulation of genes.""")
    
    submission_notes = forms.CharField(
        label="Submission notes",
        required=False, 
        help_text="""Include any additional details about the submission. For
        instance, you might indicate the approximate number of sites reported,
        whether high-throughput techniques are used or any other factor that
        might help prioritize curation.""",
        widget=forms.Textarea)
