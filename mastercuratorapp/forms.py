from django import forms
from collectfapp import curationform

class EditCurationForm(curationform.GenomeForm,
                       curationform.TechniquesForm,
                       curationform.CurationReviewForm):
    pass
