from django import forms
from collectfapp import curationform

class EditCurationForm(curationform.GenomeForm,
                       curationform.TechniquesForm,
                       curationform.CurationReviewForm):

    def clean(self):
        cleaned_data = super(EditCurationForm, self).clean()
        print 'cleaned_data', cleaned_data
        return cleaned_data

        
