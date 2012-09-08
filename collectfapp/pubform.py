from django import forms
from models import Publication
import bioutils
import sutils

def citation(pubrec):
    """Create citation string from pubmed record retrieved from NCBI"""
    title = pubrec["Title"]
    authors = pubrec["AuthorList"]
    journal = pubrec["FullJournalName"]
    return '|'.join([title, ','.join(authors), journal])

def make_pub(pubrec, cd):
    """Given publication record, retrieved from NCBI db, and entered user data,
    create models.Publication object and return it.
    """
    pmid = cd["pmid"]
    p = Publication(publication_type="pubmed",
                    pmid=pmid,
                    citation=citation(pubrec),
                    url="http://www.ncbi.nlm.nih.gov/pubmed?term=%s" % pmid,
                    pdf=None,
                    contains_promoter_data=cd["contains_promoter_data"],
                    contains_expression_data=cd["contains_expression_data"],
                    submission_notes=cd["submission_notes"],
                    curation_complete=False)
    return p
                    
class PubmedPublicationForm(forms.Form):
    pmid = forms.CharField(required=True)
    contains_promoter_data = forms.BooleanField(required=False)
    contains_expression_data = forms.BooleanField(required=False)
    submission_notes = forms.CharField(widget=forms.Textarea, required=False)

    def clean(self):         # clean fields
        cd = self.cleaned_data
        try:  # to get publication from db
            p = Publication.objects.get(pmid=cd["pmid"])
            raise forms.ValidationError("Publication is already in database.")
        except Publication.DoesNotExist:
            pubrec = bioutils.get_pubmed(cd["pmid"])
            if not pubrec:
                msg = "Cannot grab publication details from NCBI (pmid is invalid)"
                raise forms.ValidationError(msg)
            # else, add pub to database
            p = make_pub(pubrec, cd)
            sutils.sput(self.request.session, "publication", p)
        return cd
    
    

            
