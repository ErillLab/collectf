"""Forms for adding experimental techniques"""

from django import forms
import base.models as models

class AddTechniqueForm(forms.Form):
    name = forms.CharField(label="Name")
    description = forms.CharField(label="Description", widget=forms.Textarea)
    type = forms.ChoiceField(choices = models.ExperimentalTechnique.FUNCTION_CATEGORIES)
    categories = forms.ModelMultipleChoiceField(queryset=models.ExperimentalTechniqueCategory.objects.order_by('name'))

 # the foll. was added by Dinara
    EO_term = forms.CharField(label="ECO term", widget=forms.TextInput(attrs={'class':"bp_form_complete-ECO-uri", 
    																			'name':"a",'size':"100"}))