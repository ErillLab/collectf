from django import forms
import models
from django.utils.safestring import mark_safe

class BrowseForm(forms.Form):
    TF = forms.ModelChoiceField(queryset=models.TF.objects.all(),
                                label="TF",
                                initial=0)
    species = forms.ModelChoiceField(queryset=models.Strain.objects.all(),
                                     label="species",
                                     initial=0)

    browse_motif_associated_sites = forms.BooleanField(required=False)
    browse_not_motif_associated_sites = forms.BooleanField(required=False)


    
